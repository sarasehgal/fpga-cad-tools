// ece1387 a2 - analytical placer (clique model, solved w umfpack) + flow based spreading (darav fpga'19)
// to run it: ./placer <circuit> [-nogfx] [-psi lin|half|lin2|from5|dbl|big] [-weak w] [-strong w] [-seed s] [-svg] [-png] [-sweep]
// (no -psi -> tries all the psi schedules n uses the one w the smallest displacement for the anchors)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <string>
#include <vector>
#include "umfpack.h"
#ifndef NOGFX
#include "ezgl/application.hpp"
#include "ezgl/graphics.hpp"
#endif
using namespace std;

const int N = 35;  // area is 0,0 -> 35,35, unit bins so 35 x 35 = 1225 bins
typedef vector<double> Vec;

// w = width (1,2,3), all cells r 1 tall. fixed ones get fx,fy from the 2nd part of the file
struct Cell { int id, w; bool fixed = false; double fx = 0, fy = 0; };
vector<Cell> cells;
vector<vector<int>> nets;  // net -> cells on it (no repeats)
vector<int> mov, midx;     // movable cells, n cell -> its row in the matrix (-1 = fixed)

// artificial 2 pin net from a cell to a fixed point (part 3ii)
struct Anchor { int c; double x, y, w; };

// ---------- reading the file ----------
bool readfile(const char* f) {
    ifstream in(f);
    if (!in) return false;
    map<int, int> byid;          // block num -> index in cells
    map<int, vector<int>> netm;  // net num -> cells
    int b;
    while (in >> b && b != -1) {
        Cell c; c.id = b; in >> c.w;
        int n;
        while (in >> n && n != -1) netm[n].push_back(cells.size());
        byid[b] = cells.size();
        cells.push_back(c);
    }
    double x, y;
    while (in >> b && b != -1) {
        in >> x >> y;
        Cell& c = cells[byid[b]];
        c.fixed = true; c.fx = x; c.fy = y;
    }
    for (auto& kv : netm) {
        vector<int> v = kv.second;
        sort(v.begin(), v.end());
        v.erase(unique(v.begin(), v.end()), v.end());  // same block listed twice on a net = 1 pin
        if (v.size() >= 2) nets.push_back(v);
    }
    midx.assign(cells.size(), -1);
    for (int i = 0; i < (int)cells.size(); i++)
        if (!cells[i].fixed) { midx[i] = mov.size(); mov.push_back(i); }
    return true;
}

// ---------- the AP solve ----------
// X,Y hold positions for every cell. fixed ones r read from them, movable ones get overwritten
// clique model: p pin net -> p(p-1)/2 edges w weight 2/p
// for x: Q x = bx, Q = laplacian over the movable cells, fixed neighbours go to the diagonal n to b
bool solve(Vec& X, Vec& Y, const vector<Anchor>& anc) {
    int m = mov.size();
    vector<int> Ti, Tj;
    Vec Tx, bx(m, 0), by(m, 0);
    auto add = [&](int i, int j, double v) { Ti.push_back(i); Tj.push_back(j); Tx.push_back(v); };
    auto edge = [&](int a, int b, double w) {
        int ia = midx[a], ib = midx[b];
        if (ia >= 0 && ib >= 0) { add(ia, ia, w); add(ib, ib, w); add(ia, ib, -w); add(ib, ia, -w); }
        else if (ia >= 0) { add(ia, ia, w); bx[ia] += w * X[b]; by[ia] += w * Y[b]; }
        else if (ib >= 0) { add(ib, ib, w); bx[ib] += w * X[a]; by[ib] += w * Y[a]; }
        // both fixed -> nothing to solve
    };
    for (auto& net : nets) {
        double w = 2.0 / net.size();
        for (int a = 0; a < (int)net.size(); a++)
            for (int b = a + 1; b < (int)net.size(); b++) edge(net[a], net[b], w);
    }
    for (auto& an : anc) {
        int i = midx[an.c];
        add(i, i, an.w); bx[i] += an.w * an.x; by[i] += an.w * an.y;
    }
    // triplets -> compressed column (umfpack adds up the duplicates for us)
    int nz = Tx.size();
    vector<int> Ap(m + 1), Ai(nz);
    Vec Ax(nz), sx(m), sy(m);
    if (umfpack_di_triplet_to_col(m, m, nz, Ti.data(), Tj.data(), Tx.data(), Ap.data(), Ai.data(), Ax.data(), nullptr) != UMFPACK_OK) return false;
    void *sym, *num;
    umfpack_di_symbolic(m, m, Ap.data(), Ai.data(), Ax.data(), &sym, nullptr, nullptr);
    int st = umfpack_di_numeric(Ap.data(), Ai.data(), Ax.data(), sym, &num, nullptr, nullptr);
    umfpack_di_free_symbolic(&sym);
    if (st != UMFPACK_OK) { cout << "umfpack numeric status " << st << " (singular?)\n"; umfpack_di_free_numeric(&num); return false; }
    // same matrix for x n y, only b changes -> factor once solve twice
    umfpack_di_solve(UMFPACK_A, Ap.data(), Ai.data(), Ax.data(), sx.data(), bx.data(), num, nullptr, nullptr);
    umfpack_di_solve(UMFPACK_A, Ap.data(), Ai.data(), Ax.data(), sy.data(), by.data(), num, nullptr, nullptr);
    umfpack_di_free_numeric(&num);
    for (int i = 0; i < m; i++) { X[mov[i]] = sx[i]; Y[mov[i]] = sy[i]; }
    return true;
}

// start positions: fixed cells from the file, movable dont matter (solve overwrites)
void filepos(Vec& X, Vec& Y) {
    X.assign(cells.size(), 0); Y.assign(cells.size(), 0);
    for (int i = 0; i < (int)cells.size(); i++) if (cells[i].fixed) { X[i] = cells[i].fx; Y[i] = cells[i].fy; }
}

// hpwl, pins at the cell centres. only real nets (anchors r never in nets)
double hpwl(const Vec& X, const Vec& Y) {
    double t = 0;
    for (auto& net : nets) {
        double x0 = 1e18, x1 = -1e18, y0 = 1e18, y1 = -1e18;
        for (int c : net) { x0 = min(x0, X[c]); x1 = max(x1, X[c]); y0 = min(y0, Y[c]); y1 = max(y1, Y[c]); }
        t += x1 - x0 + y1 - y0;
    }
    return t;
}

// every movable cell needs a path to some fixed cell or the matrix is singular -> check first
bool connected() {
    vector<vector<int>> cn(cells.size());
    for (int k = 0; k < (int)nets.size(); k++) for (int c : nets[k]) cn[c].push_back(k);
    vector<bool> seen(cells.size(), false), ns(nets.size(), false);
    vector<int> st;
    for (int i = 0; i < (int)cells.size(); i++) if (cells[i].fixed) { seen[i] = true; st.push_back(i); }
    while (!st.empty()) {
        int c = st.back(); st.pop_back();
        for (int k : cn[c]) if (!ns[k]) { ns[k] = true; for (int d : nets[k]) if (!seen[d]) { seen[d] = true; st.push_back(d); } }
    }
    return count(seen.begin(), seen.end(), false) == 0;
}

// ---------- part 2: random io spots on the perimeter ----------
// 140 integer spots around the edge (35 per side), same kind of spots as the files use
void randomio(Vec& X, Vec& Y, mt19937& rng) {
    uniform_int_distribution<int> d(0, 4 * N - 1);
    for (int i = 0; i < (int)cells.size(); i++) {
        if (!cells[i].fixed) continue;
        int t = d(rng);
        if (t < N) { X[i] = t; Y[i] = 0; }
        else if (t < 2 * N) { X[i] = N; Y[i] = t - N; }
        else if (t < 3 * N) { X[i] = 3 * N - t; Y[i] = N; }
        else { X[i] = 0; Y[i] = 4 * N - t; }
    }
}

// ---------- part 3: flow based spreading ----------
// bins r 1x1, capacity 1. a w wide cell = w unit pieces side by side, each piece counts in the bin its centre is in.
// moves r always exactly 1 unit (up/down/left/right) so a cell keeps its offset inside the bin (quantized moves).
// mixed size trick: moving a w wide cell left/right by 1 only frees its last piece's bin n fills the bin past its
// other end -> its 1 unit of area jumping w bins, so it works like a normal edge in the path.
// moving a wide cell up/down frees w bins n fills w bins -> the path can only carry 1 unit along, so this move is
// only ok if at most 1 of the w bins it lands in is alrdy full (the path keeps going from that one, or ends if none).
Vec px, py;                  // current (spread) centres
Vec ox, oy;                  // solved centres (to measure displacement)
vector<int> use;             // pieces in each bin
vector<vector<int>> inb;     // cells w a piece in each bin
int minw = 1;                // only cells at least this wide get counted (phase 1 = wide cells only)
const int DX[4] = {1, -1, 0, 0}, DY[4] = {0, 0, 1, -1};

int clampi(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }
// bins the cell covers if its centre was at (cx,cy). bin id = row * N + col
void binsof(int c, double cx, double cy, vector<int>& out) {
    out.clear();
    int w = cells[c].w, r = clampi((int)floor(cy), 0, N - 1);
    double l = cx - w / 2.0;
    for (int k = 0; k < w; k++) out.push_back(r * N + clampi((int)floor(l + k + 0.5), 0, N - 1));
}
void place(int c, int s) {  // s = +1 add, -1 remove
    if (cells[c].w < minw) return;
    vector<int> b;
    binsof(c, px[c], py[c], b);
    for (int v : b) {
        use[v] += s;
        if (s > 0) inb[v].push_back(c);
        else inb[v].erase(find(inb[v].begin(), inb[v].end(), c));
    }
}
double disp(int c, double x, double y) { return fabs(x - ox[c]) + fabs(y - oy[c]); }
double totaldisp() { double t = 0; for (int c : mov) t += disp(c, px[c], py[c]); return t; }
int overflow() { int t = 0; for (int u : use) t += max(0, u - 1); return t; }
// stays inside 0..35?
bool inside(int c, double x, double y) {
    double h = cells[c].w / 2.0;
    return x - h > -1e-9 && x + h < N + 1e-9 && y - 0.5 > -1e-9 && y + 0.5 < N + 1e-9;
}

// one move in the path search: cell c goes 1 unit in dir d
struct Mv { int c, d, from; };

// shortest (cheapest) path from overfull bin src to a bin w room. edge cost = displacement the moved cell ends up w,
// n a move is only allowed if that displacement <= psi (the darav max movement thing).
// dijkstra over the bins, stops once the best finished path is cheaper than anything left in the heap
bool findpath(int src, double psi) {
    vector<double> dist(N * N, 1e18);
    vector<Mv> par(N * N, {-1, -1, -1});
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
    dist[src] = 0; pq.push({0, src});
    double best = 1e18; Mv bestmv = {-1, -1, -1};
    vector<int> ob, nb;
    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;
        if (d >= best) break;
        // cells alrdy moved on the path to u cant move again. soak = empty bins a wide cell going up/down on
        // the path alrdy landed in (so theyre not empty anymore for the rest of this path)
        set<int> used, soak;
        for (int v = u; v != src; v = par[v].from) {
            Mv& m = par[v];
            used.insert(m.c);
            binsof(m.c, px[m.c], py[m.c], ob);
            binsof(m.c, px[m.c] + DX[m.d], py[m.c] + DY[m.d], nb);
            for (int b : nb) if (use[b] == 0 && find(ob.begin(), ob.end(), b) == ob.end()) soak.insert(b);
        }
        for (int c : inb[u]) {
            if (used.count(c)) continue;
            binsof(c, px[c], py[c], ob);
            for (int k = 0; k < 4; k++) {
                double nx = px[c] + DX[k], ny = py[c] + DY[k];
                if (!inside(c, nx, ny)) continue;
                double cost = disp(c, nx, ny);
                if (cost > psi + 1e-9) continue;
                binsof(c, nx, ny, nb);
                // lost = old - new, gained = new - old
                vector<int> lost, gain;
                for (int v : ob) if (find(nb.begin(), nb.end(), v) == nb.end()) lost.push_back(v);
                for (int v : nb) if (find(ob.begin(), ob.end(), v) == ob.end()) gain.push_back(v);
                if (find(lost.begin(), lost.end(), u) == lost.end()) continue;  // has to actually take area out of u
                double nd = d + cost;
                // full = gained bins that alrdy hv something in them
                vector<int> full;
                bool dead = false;
                for (int v : gain) {
                    if (soak.count(v)) dead = true;  // wud double up a bin we alrdy filled on this path
                    if (use[v] >= 1) full.push_back(v);
                }
                if (dead) continue;
                if (full.empty()) {  // lands somewhere empty -> finished path
                    if (nd < best) { best = nd; bestmv = {c, k, u}; }
                } else if (full.size() == 1 && full[0] != src && nd < dist[full[0]]) {
                    // exactly 1 full bin -> the extra unit keeps going from there (for a wide cell going up/down
                    // the other bins it lands in r empty so they just soak it up)
                    dist[full[0]] = nd; par[full[0]] = {c, k, u}; pq.push({nd, full[0]});
                }
            }
        }
    }
    if (bestmv.c < 0) return false;
    // walk back n do all the moves
    vector<Mv> path = {bestmv};
    for (int v = bestmv.from; v != src; v = par[v].from) path.push_back(par[v]);
    int before = overflow();
    auto domoves = [&](int s) {
        for (auto& mv : path) { place(mv.c, -1); px[mv.c] += s * DX[mv.d]; py[mv.c] += s * DY[mv.d]; place(mv.c, +1); }
    };
    domoves(+1);
    // shudnt happen anymore but just in case a path didnt actually help -> undo it
    if (overflow() >= before) { domoves(-1); return false; }
    return true;
}

// psi schedules to try (part 3i)
struct Sched { string name; double psi0; double add, mul; };
vector<Sched> scheds = {
    {"lin", 1, 1, 1},     // 1, 2, 3, ...
    {"half", 0.5, 0.5, 1}, // 0.5, 1, 1.5, ...
    {"lin2", 2, 2, 1},    // 2, 4, 6, ...
    {"from5", 5, 1, 1},   // 5, 6, 7, ... (start big-ish)
    {"dbl", 1, 0, 2},     // 1, 2, 4, 8, ...
    {"big", 1e9, 0, 1},   // no limit at all
};

struct Snap { string name; Vec X, Y; double hp; string tag = ""; };  // tag = file name for -svg
vector<Snap> stages;  // everything the gui can show

// spreads starting from the solved spot (X,Y). leaves result in px,py. returns # outer iterations
int spread(const Vec& X, const Vec& Y, const Sched& s, bool snaps, double& lastpsi) {
    ox = X; oy = Y; px = X; py = Y;
    // first push anything hanging off the edge back in, whole units only so the offset stays the same
    for (int c : mov) {
        double h = cells[c].w / 2.0;
        while (px[c] - h < -1e-9) px[c] += 1;
        while (px[c] + h > N + 1e-9) px[c] -= 1;
        while (py[c] - 0.5 < -1e-9) py[c] += 1;
        while (py[c] + 0.5 > N + 1e-9) py[c] -= 1;
    }
    // 2 phases: first only the wide (w 2,3) cells, so they get spread out among themselves while there is still
    // space around them, then everything. (doing it all at once, the w1 cells fill up the space around a clump of
    // wide cells n then the wide ones cant go up/down anymore bc that needs w free bins at once)
    int it = 0;
    for (minw = 2; minw >= 1; minw--) {
        use.assign(N * N, 0); inb.assign(N * N, {});
        for (int c : mov) place(c, +1);
        double psi = s.psi0;
        while (overflow() > 0 && it < 500) {
            it++;
            // worst bins first
            vector<int> ov;
            for (int b = 0; b < N * N; b++) if (use[b] > 1) ov.push_back(b);
            stable_sort(ov.begin(), ov.end(), [](int a, int b) { return use[a] > use[b]; });
            int moved = 0;
            for (int b : ov)
                while (use[b] > 1 && findpath(b, psi)) moved++;
            lastpsi = psi;
            if (snaps) {
                Vec SXv = X, SYv = Y;
                for (int c : mov) { SXv[c] = px[c]; SYv[c] = py[c]; }
                char buf[150];
                snprintf(buf, sizeof buf, "spreading (%s) %s iter %d psi=%g overflow=%d", s.name.c_str(),
                         minw == 2 ? "wide cells" : "all cells", it, psi, overflow());
                stages.push_back({buf, SXv, SYv, hpwl(SXv, SYv)});
            }
            if (!moved && psi > 4 * N) break;  // psi alrdy bigger than the chip n still stuck -> give up
            psi = psi * s.mul + s.add;
        }
    }
    minw = 1;
    return it;
}

// ---------- svg dump (for the report, same look as the gui) ----------
bool shownets = true;
void svg(const Snap& s, const string& fn) {
    double k = 20;  // px per unit
    FILE* f = fopen(fn.c_str(), "w");
    if (!f) return;
    auto X = [&](double x) { return 20 + x * k; };
    auto Y = [&](double y) { return 40 + (N - y) * k; };
    fprintf(f, "<svg xmlns='http://www.w3.org/2000/svg' width='%g' height='%g' font-family='Helvetica'>\n", 40 + N * k, 75 + N * k);
    fprintf(f, "<rect width='100%%' height='100%%' fill='white'/>\n");
    fprintf(f, "<text x='20' y='25' font-size='15'>%s   HPWL=%.2f</text>\n", s.name.c_str(), s.hp);
    for (int i = 0; i <= N; i++) {
        fprintf(f, "<line x1='%g' y1='%g' x2='%g' y2='%g' stroke='#e6e6e6'/>\n", X(i), Y(0), X(i), Y(N));
        fprintf(f, "<line x1='%g' y1='%g' x2='%g' y2='%g' stroke='#e6e6e6'/>\n", X(0), Y(i), X(N), Y(i));
    }
    fprintf(f, "<rect x='%g' y='%g' width='%g' height='%g' fill='none' stroke='black' stroke-width='2'/>\n", X(0), Y(N), N * k, N * k);
    if (shownets)
        for (auto& net : nets)
            for (int a = 0; a < (int)net.size(); a++)
                for (int b = a + 1; b < (int)net.size(); b++)
                    fprintf(f, "<line x1='%.2f' y1='%.2f' x2='%.2f' y2='%.2f' stroke='rgb(90,90,90)' stroke-opacity='0.25' stroke-width='0.7'/>\n",
                            X(s.X[net[a]]), Y(s.Y[net[a]]), X(s.X[net[b]]), Y(s.Y[net[b]]));
    const char* col[] = {"", "#8ec5ff", "#7fd88f", "#ffb347"};
    for (int c = 0; c < (int)cells.size(); c++) {
        double w = cells[c].fixed ? 0.6 : cells[c].w, h = cells[c].fixed ? 0.6 : 1;
        fprintf(f, "<rect x='%.2f' y='%.2f' width='%.2f' height='%.2f' fill='%s' fill-opacity='0.85' stroke='black' stroke-width='0.6'/>\n",
                X(s.X[c] - w / 2), Y(s.Y[c] + h / 2), w * k, h * k, cells[c].fixed ? "#b22222" : col[cells[c].w]);
        fprintf(f, "<text x='%.2f' y='%.2f' font-size='8' text-anchor='middle'>%d</text>\n", X(s.X[c]), Y(s.Y[c]) + 3, cells[c].id);
    }
    // legend
    fprintf(f, "<text x='20' y='%g' font-size='12'>", 62 + N * k);
    fprintf(f, "<tspan fill='#3a8ee6'>&#9632; w=1</tspan>  <tspan fill='#2e9e45'>&#9632; w=2</tspan>  <tspan fill='#e08a00'>&#9632; w=3</tspan>  <tspan fill='#b22222'>&#9632; fixed</tspan></text>\n");
    fprintf(f, "</svg>\n");
    fclose(f);
}

// ---------- gui ----------
#ifndef NOGFX
int cur = 0;
ezgl::color wcol[] = {ezgl::BLACK, {142, 197, 255}, {127, 216, 143}, {255, 179, 71}};

void draw(ezgl::renderer* g) {
    Snap& s = stages[cur];
    // bins
    g->set_line_width(1);
    g->set_color(230, 230, 230);
    for (int i = 0; i <= N; i++) { g->draw_line({(double)i, 0}, {(double)i, (double)N}); g->draw_line({0, (double)i}, {(double)N, (double)i}); }
    g->set_color(ezgl::BLACK);
    g->set_line_width(2);
    g->draw_rectangle({0, 0}, {(double)N, (double)N});
    // rats nest = the clique edges
    if (shownets) {
        g->set_line_width(1);
        g->set_color(90, 90, 90, 70);
        for (auto& net : nets)
            for (int a = 0; a < (int)net.size(); a++)
                for (int b = a + 1; b < (int)net.size(); b++)
                    g->draw_line({s.X[net[a]], s.Y[net[a]]}, {s.X[net[b]], s.Y[net[b]]});
    }
    // cells, colour by width, fixed ones small n dark red
    g->set_font_size(9);
    for (int c = 0; c < (int)cells.size(); c++) {
        double w = cells[c].fixed ? 0.6 : cells[c].w, h = cells[c].fixed ? 0.6 : 1;
        ezgl::point2d ll(s.X[c] - w / 2, s.Y[c] - h / 2);
        if (cells[c].fixed) g->set_color(178, 34, 34, 230);
        else g->set_color(wcol[cells[c].w], 215);
        g->fill_rectangle(ll, w, h);
        g->set_color(ezgl::BLACK);
        g->set_line_width(1);
        g->draw_rectangle(ll, w, h);
        g->draw_text({s.X[c], s.Y[c]}, to_string(cells[c].id), w, h);  // only drawn if it fits (zoom in)
    }
}

void showmsg(ezgl::application* app) {
    char buf[300];
    snprintf(buf, sizeof buf, "[%d/%d] %s   HPWL = %.2f", cur + 1, (int)stages.size(), stages[cur].name.c_str(), stages[cur].hp);
    app->update_message(buf);
    app->refresh_drawing();
}
void nextcb(GtkWidget*, ezgl::application* app) { if (cur + 1 < (int)stages.size()) cur++; showmsg(app); }
void prevcb(GtkWidget*, ezgl::application* app) { if (cur > 0) cur--; showmsg(app); }
void netscb(GtkWidget*, ezgl::application* app) { shownets = !shownets; showmsg(app); }
void pngcb(GtkWidget*, ezgl::application* app) {
    string fn = "stage" + to_string(cur + 1) + ".png";
    app->get_canvas("MainCanvas")->print_png(fn.c_str());
    app->update_message("saved " + fn);
}
// -png: saves the main stages as png (thru ezgl's own print_png) n quits. runs once the window is up
string pngbase;
gboolean autosave(gpointer data) {
    auto app = (ezgl::application*)data;
    for (int i = 0; i < (int)stages.size(); i++) {
        if (stages[i].tag.empty()) continue;
        cur = i;
        string fn = pngbase + "_" + stages[i].tag + ".png";
        app->get_canvas("MainCanvas")->print_png(fn.c_str(), 900, 900);
        cout << "saved " << fn << endl;
    }
    app->quit();
    return FALSE;
}
void setup(ezgl::application* app, bool) {
    int row = 6;  // ezgl's own zoom/pan buttons use the first rows
    app->create_button("Next stage", row++, nextcb);
    app->create_button("Prev stage", row++, prevcb);
    app->create_button("Nets on/off", row++, netscb);
    app->create_button("Save PNG", row++, pngcb);
    showmsg(app);
    if (!pngbase.empty()) g_timeout_add(800, autosave, app);
}
#endif

int main(int argc, char** argv) {
    if (argc < 2) { cout << "usage: ./placer <circuit> [-nogfx] [-psi lin|half|lin2|from5|dbl|big] [-weak w] [-strong w] [-seed s] [-svg] [-png] [-sweep]\n"; return 1; }
    bool gfx = true, dosvg = false, sweep = false, dopng = false;
    setvbuf(stdout, nullptr, _IOLBF, 0);  // so the output shows up right away even when piped
    string onlypsi;
    double weak = 0.1, strong = 30;
    unsigned seed = 1387;
    for (int i = 2; i < argc; i++) {
        string a = argv[i];
        if (a == "-nogfx") gfx = false;
        else if (a == "-svg") dosvg = true;
        else if (a == "-sweep") sweep = true;
        else if (a == "-png") dopng = true;
        else if (a == "-psi" && i + 1 < argc) onlypsi = argv[++i];
        else if (a == "-weak" && i + 1 < argc) weak = atof(argv[++i]);
        else if (a == "-strong" && i + 1 < argc) strong = atof(argv[++i]);
        else if (a == "-seed" && i + 1 < argc) seed = atoi(argv[++i]);
    }
    if (!readfile(argv[1])) { cout << "cant open " << argv[1] << "\n"; return 1; }
    string name = argv[1];
    name = name.substr(name.find_last_of('/') + 1);
    int w[4] = {0, 0, 0, 0};
    for (int c : mov) w[cells[c].w]++;
    printf("%s: %d cells (%d movable: %d w1, %d w2, %d w3), %d fixed, %d nets\n", name.c_str(), (int)cells.size(),
           (int)mov.size(), w[1], w[2], w[3], (int)(cells.size() - mov.size()), (int)nets.size());
    if (!connected()) cout << "warning: some movable cells dont connect to any fixed cell, matrix will be singular\n";

    // part 1
    Vec X, Y;
    filepos(X, Y);
    if (!solve(X, Y, {})) { cout << "solve failed\n"; return 1; }
    Vec solX = X, solY = Y;
    double hp1 = hpwl(X, Y);
    printf("\n[part 1] clique model, file io spots, no spreading: HPWL = %.2f\n", hp1);
    stages.push_back({"part 1: solved (clique, no spreading)", X, Y, hp1, "p1_solved"});

    // part 2: 20 random perimeter io placements
    mt19937 rng(seed);
    Vec hps;
    for (int t = 0; t < 20; t++) {
        Vec RX = X, RY = Y;
        randomio(RX, RY, rng);
        solve(RX, RY, {});
        hps.push_back(hpwl(RX, RY));
        if (t == 0) stages.push_back({"part 2: random io trial 1", RX, RY, hps.back(), "p2_randomio1"});
    }
    double mean = 0, sd = 0;
    for (double h : hps) mean += h;
    mean /= hps.size();
    for (double h : hps) sd += (h - mean) * (h - mean);
    sd = sqrt(sd / (hps.size() - 1));  // sample std dev
    printf("\n[part 2] 20 random io placements (seed %u):\n  ", seed);
    for (double h : hps) printf("%.1f ", h);
    printf("\n  mean HPWL = %.2f, std dev = %.2f (min %.2f, max %.2f)\n", mean, sd, *min_element(hps.begin(), hps.end()),
           *max_element(hps.begin(), hps.end()));

    // part 3i: try the psi schedules
    printf("\n[part 3i] flow based spreading, starting from the part 1 solve\n");
    printf("  %-6s %6s %8s %10s %9s %9s %9s\n", "psi", "iters", "lastpsi", "totaldisp", "maxdisp", "HPWL", "overflow");
    string bestname;
    double bestd = 1e18;
    Sched best;
    for (auto& s : scheds) {
        if (!onlypsi.empty() && s.name != onlypsi) continue;
        double lp = 0;
        int it = spread(solX, solY, s, false, lp);
        Vec SXv = solX, SYv = solY;
        double md = 0;
        for (int c : mov) { SXv[c] = px[c]; SYv[c] = py[c]; md = max(md, disp(c, px[c], py[c])); }
        double td = totaldisp();
        printf("  %-6s %6d %8g %10.2f %9.2f %9.2f %9d\n", s.name.c_str(), it, lp > 1e8 ? -1 : lp, td, md, hpwl(SXv, SYv), overflow());
        if (overflow() == 0 && td < bestd) { bestd = td; best = s; }
    }
    printf("  best (smallest total displacement): %s\n", best.name.c_str());

    // redo the best one w snapshots for the gui
    double lp;
    spread(solX, solY, best, true, lp);
    Vec SXv = solX, SYv = solY;
    for (int c : mov) { SXv[c] = px[c]; SYv[c] = py[c]; }
    double hps3 = hpwl(SXv, SYv);
    char b3[150];
    snprintf(b3, sizeof b3, "part 3i: spread (psi %s), total disp %.0f", best.name.c_str(), totaldisp());
    stages.push_back({b3, SXv, SYv, hps3, "p3i_spread"});

    // part 3ii: anchors at the spread spots, weak n strong
    printf("\n[part 3ii] anchors at the spread spots (psi schedule %s)\n", best.name.c_str());
    printf("  solved HPWL %.2f -> spread HPWL %.2f\n", hp1, hps3);
    for (double aw : {weak, strong}) {
        vector<Anchor> anc;
        for (int c : mov) anc.push_back({c, px[c], py[c], aw});
        Vec AX = solX, AY = solY;
        solve(AX, AY, anc);
        double d = 0;
        for (int c : mov) d += fabs(AX[c] - px[c]) + fabs(AY[c] - py[c]);
        double h = hpwl(AX, AY);
        printf("  anchor weight %-6g -> HPWL %.2f  (avg dist from anchor %.3f)\n", aw, h, d / mov.size());
        char buf[100];
        snprintf(buf, sizeof buf, "part 3ii: re-solved w %s anchors (w=%g)", aw == weak ? "weak" : "strong", aw);
        stages.push_back({buf, AX, AY, h, aw == weak ? "p3ii_weak" : "p3ii_strong"});
    }

    // extra: a sweep of anchor weights to see the trend (not needed for the handout, just for the report)
    if (sweep) {
        printf("  anchor weight sweep:\n");
        for (double aw : {0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0}) {
            vector<Anchor> anc;
            for (int c : mov) anc.push_back({c, px[c], py[c], aw});
            Vec AX = solX, AY = solY;
            solve(AX, AY, anc);
            double d = 0;
            for (int c : mov) d += fabs(AX[c] - px[c]) + fabs(AY[c] - py[c]);
            printf("    w=%-6g HPWL %9.2f  total dist from anchors %8.2f\n", aw, hpwl(AX, AY), d);
        }
    }

    if (dosvg) {
        string base = name.substr(0, name.find('.'));
        for (auto& s : stages)
            if (!s.tag.empty()) svg(s, base + "_" + s.tag + ".svg");
    }

    if (!gfx) return 0;
#ifndef NOGFX
    ezgl::application::settings st;
    st.main_ui_resource = "main.ui";
    st.window_identifier = "MainWindow";
    st.canvas_identifier = "MainCanvas";
    if (dopng) pngbase = name.substr(0, name.find('.'));
    ezgl::application app(st);
    app.add_canvas("MainCanvas", draw, ezgl::rectangle({-1, -1}, N + 2, N + 2));
    return app.run(setup, nullptr, nullptr, nullptr);
#else
    (void)dopng;
    return 0;
#endif
}
