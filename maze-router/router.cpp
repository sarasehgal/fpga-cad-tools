// ece1387 a1 - fpga maze router (lee-moore)
// to run it: ./router <circuit> [-fc] [-fs4] [-minw] [-nogfx]

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <vector>
#include "ezgl/application.hpp"
#include "ezgl/graphics.hpp"
using namespace std;

// pin = which block + which of the 8 pins
struct Pin { int x, y, p; };
// one src can drive multiple sinks -> keep them together so they can share wire
struct Net { Pin src; vector<Pin> sinks; };

int n, W;                         // n x n blocks, W tracks per channel
bool fcodd = false, fs4 = false;  // part 4 n part 5
vector<Net> nets;
vector<vector<int>> adj;          // seg -> segs it can switch to
vector<int> owner;                // which net uses each seg, -1 = free
vector<pair<Pin, int>> stubs;     // pin -> seg its hooked to (for drawing)
vector<pair<int, int>> sw;        // seg pairs joined by a switch (for drawing)

// seg ids - horiz ones first then vert. hseg(x,y) is btwn sb(x,y) n sb(x+1,y), vseg(x,y) btwn sb(x,y) n sb(x,y+1)
// n+1 horiz channels (one under each row + top one), each has n segs, each seg has W tracks
int nh() { return (n + 1) * n * W; }
int hseg(int x, int y, int t) { return (y * n + x) * W + t; }
int vseg(int x, int y, int t) { return nh() + (x * n + y) * W + t; }
// opposite of hseg/vseg, id -> (vert?, x, y, track). only rly need it for drawing
void info(int s, bool& v, int& x, int& y, int& t) {
    v = s >= nh();
    if (v) s -= nh();
    t = s % W; s /= W;
    if (v) { x = s / n; y = s % n; }
    else { x = s % n; y = s / n; }
}

// switches go both ways so add both directions. -1 = seg doesnt exist (edge of chip)
void link(int a, int b) {
    if (a < 0 || b < 0) return;
    adj[a].push_back(b);
    adj[b].push_back(a);
}

// builds the switch blocks
void build() {
    adj.assign(2 * nh(), {});
    // every sb is at a corner (i,j), i n j go 0..n
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++)
            for (int t = 0; t < W; t++) {
                // the 4 segs touching this sb on track t
                int lf = i > 0 ? hseg(i - 1, j, t) : -1, rt = i < n ? hseg(i, j, t) : -1;
                int dn = j > 0 ? vseg(i, j - 1, t) : -1, up = j < n ? vseg(i, j, t) : -1;
                // planar: track t only to track t on other 3 sides
                link(lf, rt); link(lf, dn); link(lf, up); link(rt, dn); link(rt, up); link(dn, up);
                // fs4: extra switch, left t <-> up t+1 n right t <-> down t+1
                // (w W=1 it wud just be t to t again so skip)
                if (fs4 && W > 1) {
                    int nt = (t + 1) % W;
                    link(lf, j < n ? vseg(i, j, nt) : -1);
                    link(rt, j > 0 ? vseg(i, j - 1, nt) : -1);
                }
            }
}

// segs a pin can reach thru its connection block
// pins 1,2 bottom, 3,4 left, 5,6 top, 7,8 right
vector<int> pinsegs(Pin p) {
    vector<int> r;
    for (int t = 0; t < W; t++) {
        // part 4 - inputs only get tracks w same odd/even as the pin. 4 n 7 r outputs so they get all
        if (fcodd && p.p != 4 && p.p != 7 && t % 2 != p.p % 2) continue;
        if (p.p <= 2) r.push_back(hseg(p.x, p.y, t));
        else if (p.p <= 4) r.push_back(vseg(p.x, p.y, t));
        else if (p.p <= 6) r.push_back(hseg(p.x, p.y + 1, t));   // top = channel above so y+1
        else r.push_back(vseg(p.x + 1, p.y, t));                  // right = channel to the right so x+1
    }
    return r;
}

// bfs from src pin + all segs net id alrdy uses, til we hit a seg next to the sink
bool routesink(int id, Pin src, Pin snk) {
    int tot = adj.size();
    vector<int> lab(tot, -1), prev(tot, -1);   // lab = ripple number, prev = who labelled me
    vector<bool> tgt(tot, false);
    // can land on a sink seg if its free or alrdy ours
    for (int s : pinsegs(snk)) if (owner[s] < 0 || owner[s] == id) tgt[s] = true;

    queue<int> q;
    // start = free segs at the src pin + everything this net alrdy has (so branches share wire)
    for (int s : pinsegs(src)) if (owner[s] < 0) { lab[s] = 0; q.push(s); }
    for (int s = 0; s < tot; s++) if (owner[s] == id) { lab[s] = 0; q.push(s); }

    while (!q.empty()) {
        int j = q.front(); q.pop();
        // checking here not when labelling, still shortest bc fifo
        if (tgt[j]) {
            // found it, trace back n mark segs used
            stubs.push_back({snk, j});
            int s = j;
            while (prev[s] >= 0) {
                owner[s] = id;
                sw.push_back({s, prev[s]});
                s = prev[s];
            }
            if (owner[s] != id) stubs.push_back({src, s});  // new branch off src pin
            owner[s] = id;
            return true;
        }
        // ripple out. skip anything labelled alrdy or used by smone
        for (int k : adj[j])
            if (lab[k] < 0 && owner[k] < 0) { lab[k] = lab[j] + 1; prev[k] = j; q.push(k); }
    }
    return false;   // ran out of places to go, cant route this one
}

// routes all nets in given order, returns net that failed or -1 if all good
int routeall(vector<int>& order) {
    owner.assign(adj.size(), -1);   // wipe everything first
    stubs.clear(); sw.clear();
    for (int i : order)
        for (Pin s : nets[i].sinks)
            if (!routesink(i, nets[i].src, s)) return i;
    return -1;
}

// if a net fails, move it to front n retry
// idea: the ones that fail r the hard ones so let them pick first
bool tryorder(vector<int> order) {
    for (int it = 0; it < 100; it++) {   // 100 is j so it doesnt loop forever
        int f = routeall(order);
        if (f < 0) return true;
        order.erase(find(order.begin(), order.end(), f));
        order.insert(order.begin(), f);
    }
    return false;
}

// try w tracks w file order, if that fails try long nets first
// (long first alone was worse on some ccts so its only a backup)
bool tryw(int w) {
    W = w;
    build();
    vector<int> order;
    for (int i = 0; i < (int)nets.size(); i++) order.push_back(i);
    if (tryorder(order)) return true;
    // len = manhattan dist to the farthest sink
    auto len = [](Net& t) { int m = 0; for (Pin s : t.sinks) m = max(m, abs(s.x - t.src.x) + abs(s.y - t.src.y)); return m; };
    stable_sort(order.begin(), order.end(), [&](int a, int b) { return len(nets[a]) > len(nets[b]); });
    return tryorder(order);
}

// segs used = anything w an owner
int usedsegs() { return count_if(owner.begin(), owner.end(), [](int o) { return o >= 0; }); }

// drawing stuff
double chw, blk, step;  // chan width, block size, chan+block

// sb(i,j) = square at (i*step, j*step) w size chw. gives end of seg s at that sb
// track t sits t+1 in from the edge of its channel
ezgl::point2d endat(int s, int i, int j) {
    bool v; int x, y, t;
    info(s, v, x, y, t);
    if (v) return {x * step + t + 1, j == y ? j * step + chw : j * step};
    return {i == x ? i * step + chw : i * step, y * step + t + 1};
}

// where to put the dot for each pin on the block edge (index 0 not used)
ezgl::point2d pinpos(Pin p) {
    double bx = p.x * step + chw, by = p.y * step + chw, a = blk / 3, b = 2 * blk / 3;
    double xs[] = {0, bx + b, bx + a, bx, bx, bx + b, bx + a, bx + blk, bx + blk};
    double ys[] = {0, by, by, by + a, by + b, by + blk, by + blk, by + b, by + a};
    return {xs[p.p], ys[p.p]};
}

// just cycle thru some colours, repeats after 11 but fine
ezgl::color netcolor(int id) {
    ezgl::color cs[] = {ezgl::RED, ezgl::BLUE, ezgl::DARK_GREEN, ezgl::ORANGE, ezgl::PURPLE, ezgl::MAGENTA,
                        ezgl::SADDLE_BROWN, ezgl::CORAL, ezgl::TURQUOISE, ezgl::FIRE_BRICK, ezgl::DARK_SLATE_BLUE};
    return cs[id % 11];
}

// whole seg, sb to sb
void drawseg(ezgl::renderer* g, int s) {
    bool v; int x, y, t;
    info(s, v, x, y, t);
    g->draw_line(endat(s, x, y), endat(s, x + !v, y + v));
}

// ezgl calls this whenever it redraws (zoom etc)
void draw(ezgl::renderer* g) {
    // blocks
    g->set_line_width(1);
    for (int x = 0; x < n; x++)
        for (int y = 0; y < n; y++) {
            g->set_color(ezgl::GREY_75);
            g->fill_rectangle({x * step + chw, y * step + chw}, blk, blk);
            g->set_color(ezgl::BLACK);
            g->draw_text({x * step + chw + blk / 2, y * step + chw + blk / 2}, to_string(x) + "," + to_string(y));
        }
    // switch blocks
    g->set_color(210, 225, 250);
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++) g->fill_rectangle({i * step, j * step}, chw, chw);
    // all tracks faint
    g->set_color(220, 220, 220);
    for (int s = 0; s < (int)adj.size(); s++) drawseg(g, s);
    // used segs
    g->set_line_width(3);
    for (int s = 0; s < (int)adj.size(); s++)
        if (owner[s] >= 0) { g->set_color(netcolor(owner[s])); drawseg(g, s); }
    // switches used (line inside the sb)
    for (auto& e : sw) {
        bool va, vb; int xa, ya, ta, xb, yb, tb;
        info(e.first, va, xa, ya, ta);
        info(e.second, vb, xb, yb, tb);
        g->set_color(netcolor(owner[e.first]));
        // each seg has 2 end sbs, find the one both segs share
        for (int p = 0; p < 2; p++)
            for (int q = 0; q < 2; q++) {
                int i = xa + p * !va, j = ya + p * va;
                if (i == xb + q * !vb && j == yb + q * vb) g->draw_line(endat(e.first, i, j), endat(e.second, i, j));
            }
    }
    // pin stubs
    g->set_line_width(2);
    for (auto& sp : stubs) {
        bool v; int x, y, t;
        info(sp.second, v, x, y, t);
        // line from the pin straight across to its track
        ezgl::point2d a = pinpos(sp.first), b = a;
        if (v) b.x = x * step + t + 1;
        else b.y = y * step + t + 1;
        g->set_color(netcolor(owner[sp.second]));
        g->draw_line(a, b);
        g->fill_arc(a, 0.25, 0, 360);
    }
}

string msg;   // W + segs, shows in the status bar
void setup(ezgl::application* app, bool) { app->update_message(msg); }

int main(int argc, char** argv) {
    if (argc < 2) { cout << "usage: ./router <circuit> [-fc] [-fs4] [-minw] [-nogfx]\n"; return 1; }
    bool minw = false, gfx = true;
    for (int i = 2; i < argc; i++) {
        string a = argv[i];
        if (a == "-fc") fcodd = true;
        if (a == "-fs4") fs4 = true;
        if (a == "-minw") minw = true;
        if (a == "-nogfx") gfx = false;   // for quick runs
    }

    // read file, group sinks by src pin into nets
    ifstream in(argv[1]);
    if (!in) { cout << "cant open " << argv[1] << "\n"; return 1; }
    int w;
    in >> n >> w;
    Pin a, b;
    // stops at the -1 -1 -1 -1 -1 -1 line
    while (in >> a.x >> a.y >> a.p >> b.x >> b.y >> b.p && a.x != -1) {
        // alrdy hv a net w this src? if not make one
        int k = 0;
        while (k < (int)nets.size() && !(nets[k].src.x == a.x && nets[k].src.y == a.y && nets[k].src.p == a.p)) k++;
        if (k == (int)nets.size()) nets.push_back({a, {}});
        nets[k].sinks.push_back(b);
    }

    bool ok;
    // minw: go up from 1 til it works (100 cap j in case)
    if (minw) {
        for (w = 1; !(ok = tryw(w)) && w < 100; w++);
    } else ok = tryw(w);

    msg = string(argv[1]) + (fcodd ? " fc=odd/even" : "") + (fs4 ? " fs=4" : "") + "  W=" + to_string(W) +
          (ok ? "  routed, segs used=" : "  FAILED, segs used=") + to_string(usedsegs());
    cout << msg << "\n";
    if (!gfx) return 0;

    // sizes for drawing, channel gets 1 unit per track
    chw = W + 1;
    blk = max(chw, 6.0);
    step = blk + chw;
    // same setup as the ezgl example
    ezgl::application::settings s;
    s.main_ui_resource = "main.ui";
    s.window_identifier = "MainWindow";
    s.canvas_identifier = "MainCanvas";
    ezgl::application app(s);
    app.add_canvas("MainCanvas", draw, ezgl::rectangle({0, 0}, n * step + chw, n * step + chw));
    return app.run(setup, nullptr, nullptr, nullptr);
}
