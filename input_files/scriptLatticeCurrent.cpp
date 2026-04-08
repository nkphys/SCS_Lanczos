#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

struct SetData {
    int set_id;
    int n;
    int m;
    bool has_sites;
    bool has_total;
    double total_real;
    double total_imag;

    SetData() {
        set_id = -1;
        n = -1;
        m = -1;
        has_sites = false;
        has_total = false;
        total_real = 0.0;
        total_imag = 0.0;
    }
};

struct BondCurrent {
    int set_id;
    int n;
    int m;
    double current;
};

struct Point {
    double x;
    double y;
};

enum TotalMode {
    TOTAL_MODE_TOTAL,
    TOTAL_MODE_QUANTUM
};

static bool ParseRunOut(const string &filepath, vector<BondCurrent> &bonds, int &max_site, TotalMode mode) {
    ifstream infile(filepath.c_str());
    if (!infile.is_open()) {
        cerr << "Unable to open input file: " << filepath << endl;
        return false;
    }

    regex set_re("^\\s*Set\\s*=\\s*([0-9]+)");
    regex sites_re("sites=\\(([-]?[0-9]+),([-]?[0-9]+),([-]?[0-9]+),([-]?[0-9]+)\\)");
    regex total_re("^\\s*Total\\s+for\\s+set\\s+([0-9]+)\\s*=\\s*\\(([-+0-9.eE]+),([-+0-9.eE]+)\\)");
    regex total_quantum_re("^\\s*Total\\s+quantum\\s+for\\s+set\\s+([0-9]+)\\s*=\\s*\\(([-+0-9.eE]+),([-+0-9.eE]+)\\)");

    map<int, SetData> set_map;
    int current_set = -1;
    string line;

    while (getline(infile, line)) {
        smatch match;

        if (regex_search(line, match, set_re)) {
            current_set = atoi(match[1].str().c_str());
            if (set_map.find(current_set) == set_map.end()) {
                SetData data;
                data.set_id = current_set;
                set_map[current_set] = data;
            }
            continue;
        }

        if (current_set >= 0 && regex_search(line, match, sites_re)) {
            if (!set_map[current_set].has_sites) {
                int n = atoi(match[1].str().c_str());
                int m = atoi(match[2].str().c_str());
                set_map[current_set].n = n;
                set_map[current_set].m = m;
                set_map[current_set].has_sites = true;
            }
            continue;
        }

        bool total_match = false;
        if (mode == TOTAL_MODE_TOTAL && regex_search(line, match, total_re)) {
            total_match = true;
        }
        if (mode == TOTAL_MODE_QUANTUM && regex_search(line, match, total_quantum_re)) {
            total_match = true;
        }

        if (total_match) {
            int set_id = atoi(match[1].str().c_str());
            double total_real = atof(match[2].str().c_str());
            double total_imag = atof(match[3].str().c_str());

            if (set_map.find(set_id) == set_map.end()) {
                SetData data;
                data.set_id = set_id;
                set_map[set_id] = data;
            }
            set_map[set_id].total_real = total_real;
            set_map[set_id].total_imag = total_imag;
            set_map[set_id].has_total = true;
        }
    }

    bonds.clear();
    max_site = -1;

    for (map<int, SetData>::const_iterator it = set_map.begin(); it != set_map.end(); ++it) {
        const SetData &data = it->second;
        if (data.has_sites && data.has_total) {
            BondCurrent item;
            item.set_id = data.set_id;
            item.n = data.n;
            item.m = data.m;
            item.current = data.total_real;
            bonds.push_back(item);
            max_site = max(max_site, max(item.n, item.m));
        }
    }

    sort(bonds.begin(), bonds.end(), [](const BondCurrent &a, const BondCurrent &b) {
        return a.set_id < b.set_id;
    });

    return true;
}

static Point UnitCellOffset(int local_index, int sites_per_unit_cell, double unit) {
    Point p;

    if (sites_per_unit_cell == 2) {
        const double off = 0.22 * unit;
        const double low = off;
        const double high = unit - off;
        if (local_index == 0) { p.x = low;  p.y = low; }
        else { p.x = low; p.y = high; }
        return p;
    }

    if (sites_per_unit_cell == 4) {
        const double off = 0.22 * unit;
        const double low = off;
        const double high = unit - off;
        if (local_index == 0) { p.x = low;  p.y = low; }
        else if (local_index == 1) { p.x = high; p.y = low; }
        else if (local_index == 2) { p.x = high; p.y = high; }
        else { p.x = low; p.y = high; }
        return p;
    }

    if (sites_per_unit_cell == 3) {
        const double off = 0.22 * unit;
        const double low = off;
        const double high = unit - off;
        const double mid = 0.5 * unit;
        if (local_index == 0) { p.x = low;  p.y = low; }
        else if (local_index == 1) { p.x = high; p.y = low; }
        else { p.x = mid; p.y = high; }
        return p;
    }

    double angle = 2.0 * M_PI * (static_cast<double>(local_index) / static_cast<double>(sites_per_unit_cell));
    double r = 0.3 * unit;
    p.x = 0.5 * unit + r * cos(angle);
    p.y = 0.5 * unit + r * sin(angle);
    return p;
}

static Point SiteToPoint(int site, int Lx, int Ly, int sites_per_unit_cell, double unit) {
    if (sites_per_unit_cell == 2) {
        int a = site % 2;
        int cell_linear = site / 2;
        int ux = cell_linear % Lx;
        int uy = cell_linear / Lx;

        double x_model = ux + ((a == 0) ? 0.5 : 0.0);
        double y_model = uy + ((a == 0) ? 0.0 : 0.5);

        Point p;
        p.x = x_model * unit;
        p.y = (Ly - y_model) * unit;
        return p;
    }

    int cell_index = site / sites_per_unit_cell;
    int local_index = site % sites_per_unit_cell;

    int x = cell_index % Lx;
    int y = cell_index / Lx;
    int y_screen = (Ly - 1 - y);

    Point local = UnitCellOffset(local_index, sites_per_unit_cell, unit);

    Point p;
    p.x = x * unit + local.x;
    p.y = y_screen * unit + local.y;
    return p;
}

static Point ShortestPeriodicDelta(const Point &from, const Point &to, double width, double height) {
    Point delta;
    delta.x = to.x - from.x;
    delta.y = to.y - from.y;

    if (delta.x > 0.5 * width) { delta.x -= width; }
    if (delta.x < -0.5 * width) { delta.x += width; }
    if (delta.y > 0.5 * height) { delta.y -= height; }
    if (delta.y < -0.5 * height) { delta.y += height; }

    return delta;
}

static pair<int, int> BestImageShift(const Point &from,
                                     const Point &to_center,
                                     double lattice_w,
                                     double lattice_h,
                                     int cushion_cells) {
    pair<int, int> best_shift = make_pair(0, 0);
    double best_d2 = 1.0e300;

    for (int sx = -cushion_cells; sx <= cushion_cells; sx++) {
        for (int sy = -cushion_cells; sy <= cushion_cells; sy++) {
            double tx = to_center.x + sx * lattice_w;
            double ty = to_center.y + sy * lattice_h;
            double dx = tx - from.x;
            double dy = ty - from.y;
            double d2 = dx * dx + dy * dy;
            if (d2 < best_d2) {
                best_d2 = d2;
                best_shift = make_pair(sx, sy);
            }
        }
    }

    return best_shift;
}

static bool WriteSvg(const string &outfile,
                     const vector<BondCurrent> &bonds,
                     int total_sites,
                     int Lx,
                     int Ly,
                     int sites_per_unit_cell,
                     double current_threshold,
                     const string &total_label) {

    ofstream svg(outfile.c_str());
    if (!svg.is_open()) {
        cerr << "Unable to open output file: " << outfile << endl;
        return false;
    }

    const double unit = 140.0;
    const double margin = 70.0;
    const int cushion_cells = 1;
    const int grid_cells_x = (2 * cushion_cells + 1) * Lx;
    const int grid_cells_y = (2 * cushion_cells + 1) * Ly;
    const double lattice_w = Lx * unit;
    const double lattice_h = Ly * unit;
    const double legend_panel_w = 260.0;
    const double extended_w = (2 * cushion_cells + 1) * lattice_w;
    const double extended_h = (2 * cushion_cells + 1) * lattice_h;
    const double canvas_w = extended_w + 2.0 * margin + legend_panel_w;
    const double canvas_h = extended_h + 2.0 * margin;
    const double origin_x = margin + cushion_cells * lattice_w;
    const double origin_y = margin + cushion_cells * lattice_h;

    vector<Point> site_pos(total_sites);
    for (int s = 0; s < total_sites; s++) {
        Point p = SiteToPoint(s, Lx, Ly, sites_per_unit_cell, unit);
        p.x += origin_x;
        p.y += origin_y;
        site_pos[s] = p;
    }

    double max_abs_current = 0.0;
    for (size_t i = 0; i < bonds.size(); i++) {
        max_abs_current = max(max_abs_current, fabs(bonds[i].current));
    }
    if (max_abs_current == 0.0) {
        max_abs_current = 1.0;
    }

    set<pair<int, int> > unique_bonds;
    for (size_t i = 0; i < bonds.size(); i++) {
        int a = min(bonds[i].n, bonds[i].m);
        int b = max(bonds[i].n, bonds[i].m);
        unique_bonds.insert(make_pair(a, b));
    }

    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << canvas_w
        << "\" height=\"" << canvas_h << "\" viewBox=\"0 0 " << canvas_w << " " << canvas_h << "\">\n";

    svg << "<defs>\n";
    svg << "<marker id=\"arrowPos\" viewBox=\"0 0 10 10\" refX=\"9\" refY=\"5\" markerWidth=\"7\" markerHeight=\"7\" orient=\"auto-start-reverse\">\n";
    svg << "<path d=\"M 0 0 L 10 5 L 0 10 z\" fill=\"#1d4ed8\"/>\n";
    svg << "</marker>\n";
    svg << "<marker id=\"arrowNeg\" viewBox=\"0 0 10 10\" refX=\"9\" refY=\"5\" markerWidth=\"7\" markerHeight=\"7\" orient=\"auto-start-reverse\">\n";
    svg << "<path d=\"M 0 0 L 10 5 L 0 10 z\" fill=\"#1d4ed8\"/>\n";
    svg << "</marker>\n";
    svg << "</defs>\n";

    svg << "<rect x=\"0\" y=\"0\" width=\"" << canvas_w << "\" height=\"" << canvas_h << "\" fill=\"white\"/>\n";

    svg << "<text x=\"" << 0.5 * canvas_w << "\" y=\"34\" text-anchor=\"middle\" font-size=\"22\" font-family=\"Arial\">Lattice currents from "
        << total_label << " i</text>\n";

    svg << "<rect x=\"" << margin << "\" y=\"" << margin << "\" width=\"" << extended_w << "\" height=\"" << extended_h
        << "\" fill=\"#fafafa\" stroke=\"#e5e7eb\" stroke-width=\"1.0\"/>\n";

    for (int gx = 0; gx <= grid_cells_x; gx++) {
        double x = margin + gx * unit;
        svg << "<line x1=\"" << x << "\" y1=\"" << margin
            << "\" x2=\"" << x << "\" y2=\"" << (margin + extended_h)
            << "\" stroke=\"#eef2f7\" stroke-width=\"1\"/>\n";
    }
    for (int gy = 0; gy <= grid_cells_y; gy++) {
        double y = margin + gy * unit;
        svg << "<line x1=\"" << margin << "\" y1=\"" << y
            << "\" x2=\"" << (margin + extended_w) << "\" y2=\"" << y
            << "\" stroke=\"#eef2f7\" stroke-width=\"1\"/>\n";
    }

    svg << "<rect x=\"" << origin_x << "\" y=\"" << origin_y << "\" width=\"" << lattice_w << "\" height=\"" << lattice_h
        << "\" fill=\"#fafafa\" stroke=\"#d1d5db\" stroke-width=\"1.2\"/>\n";

    for (set<pair<int, int> >::const_iterator it = unique_bonds.begin(); it != unique_bonds.end(); ++it) {
        int a = it->first;
        int b = it->second;
        Point pa = site_pos[a];
        Point pb_center = site_pos[b];
        pair<int, int> shift = BestImageShift(pa, pb_center, lattice_w, lattice_h, cushion_cells);
        Point pb;
        pb.x = pb_center.x + shift.first * lattice_w;
        pb.y = pb_center.y + shift.second * lattice_h;

        svg << "<line x1=\"" << pa.x << "\" y1=\"" << pa.y
            << "\" x2=\"" << pb.x << "\" y2=\"" << pb.y
            << "\" stroke=\"#cbd5e1\" stroke-width=\"1.5\"/>\n";
    }

    for (size_t i = 0; i < bonds.size(); i++) {
        const BondCurrent &bc = bonds[i];
        if (fabs(bc.current) < current_threshold) {
            continue;
        }

        Point p_n = site_pos[bc.n];
        Point p_m_center = site_pos[bc.m];
        pair<int, int> shift = BestImageShift(p_n, p_m_center, lattice_w, lattice_h, cushion_cells);
        Point p_m;
        p_m.x = p_m_center.x + shift.first * lattice_w;
        p_m.y = p_m_center.y + shift.second * lattice_h;

        Point d;
        d.x = p_m.x - p_n.x;
        d.y = p_m.y - p_n.y;

        double len = sqrt(d.x * d.x + d.y * d.y);
        if (len < 1e-12) {
            continue;
        }

        double ux = d.x / len;
        double uy = d.y / len;
        double edge_pad = 13.0;

        Point p_n_in;
        p_n_in.x = p_n.x + edge_pad * ux;
        p_n_in.y = p_n.y + edge_pad * uy;

        Point p_m_in;
        p_m_in.x = p_n.x + (len - edge_pad) * ux;
        p_m_in.y = p_n.y + (len - edge_pad) * uy;

        Point p1;
        Point p2;

        if (bc.current > 0.0) {
            p1 = p_m_in;
            p2 = p_n_in;
        }
        else {
            p1 = p_n_in;
            p2 = p_m_in;
        }

        double scaled = fabs(bc.current) / max_abs_current;
        double stroke_w = 1.5 + 7.0 * scaled;

        string color = "#1d4ed8";
        string marker = (bc.current > 0.0) ? "url(#arrowPos)" : "url(#arrowNeg)";

        int from = (bc.current > 0.0) ? bc.m : bc.n;
        int to = (bc.current > 0.0) ? bc.n : bc.m;

        for (int sx = -cushion_cells; sx <= cushion_cells; sx++) {
            for (int sy = -cushion_cells; sy <= cushion_cells; sy++) {
                Point p1_img;
                p1_img.x = p1.x + sx * lattice_w;
                p1_img.y = p1.y + sy * lattice_h;

                Point p2_img;
                p2_img.x = p2.x + sx * lattice_w;
                p2_img.y = p2.y + sy * lattice_h;

                svg << "<line x1=\"" << p1_img.x << "\" y1=\"" << p1_img.y
                    << "\" x2=\"" << p2_img.x << "\" y2=\"" << p2_img.y
                    << "\" stroke=\"" << color << "\" stroke-width=\"" << stroke_w
                    << "\" marker-end=\"" << marker << "\">\n";

                svg << "<title>set=" << bc.set_id << " n=" << bc.n << " m=" << bc.m << " total="
                    << scientific << setprecision(6) << bc.current
                    << " direction=" << from << "->" << to << "</title>\n";

                svg << "</line>\n";
            }
        }
    }

    for (int s = 0; s < total_sites; s++) {
        for (int sx = -cushion_cells; sx <= cushion_cells; sx++) {
            for (int sy = -cushion_cells; sy <= cushion_cells; sy++) {
                Point p;
                p.x = site_pos[s].x + sx * lattice_w;
                p.y = site_pos[s].y + sy * lattice_h;

                bool is_center = (sx == 0 && sy == 0);
                double r = is_center ? 7.5 : 4.8;
                string fill = is_center ? "#111827" : "#94a3b8";
                string text_fill = is_center ? "#111827" : "#64748b";
                int font_size = is_center ? 14 : 12;

                svg << "<circle cx=\"" << p.x << "\" cy=\"" << p.y
                    << "\" r=\"" << r << "\" fill=\"" << fill << "\"/>\n";
                svg << "<text x=\"" << (p.x + 8.0) << "\" y=\"" << (p.y - 7.0)
                    << "\" font-size=\"" << font_size << "\" font-family=\"Arial\" font-weight=\"700\" fill=\"" << text_fill << "\">" << s << "</text>\n";
            }
        }
    }

    double lx = margin + lattice_w + 20.0;
    double ly = margin + 10.0;
    svg << "<rect x=\"" << lx << "\" y=\"" << ly << "\" width=\"220\" height=\"88\" fill=\"#ffffff\" stroke=\"#d1d5db\"/>\n";
    svg << "<line x1=\"" << (lx + 16) << "\" y1=\"" << (ly + 20) << "\" x2=\"" << (lx + 84) << "\" y2=\"" << (ly + 20)
        << "\" stroke=\"#1d4ed8\" stroke-width=\"3\" marker-end=\"url(#arrowPos)\"/>\n";
    svg << "<text x=\"" << (lx + 94) << "\" y=\"" << (ly + 25) << "\" font-size=\"13\" font-family=\"Arial\">Total &gt; 0 : m -> n</text>\n";
    svg << "<line x1=\"" << (lx + 16) << "\" y1=\"" << (ly + 44) << "\" x2=\"" << (lx + 84) << "\" y2=\"" << (ly + 44)
        << "\" stroke=\"#1d4ed8\" stroke-width=\"3\" marker-end=\"url(#arrowNeg)\"/>\n";
    svg << "<text x=\"" << (lx + 94) << "\" y=\"" << (ly + 49) << "\" font-size=\"13\" font-family=\"Arial\">Total &lt; 0 : n -> m</text>\n";

    svg << "</svg>\n";
    svg.close();
    return true;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        cout << "Usage: ./scriptLatticeCurrent Lx Ly [sites_per_unit_cell] [Run_out.txt] [output.svg] [threshold] [total_mode]\n";
        cout << "total_mode: total | quantum (default: total)\n";
        cout << "Example: ./scriptLatticeCurrent 2 2 3 Run_out.txt lattice_current.svg 1e-12 quantum\n";
        return 1;
    }

    int Lx = atoi(argv[1]);
    int Ly = atoi(argv[2]);
    int sites_per_unit_cell = 2;
    string input_file = "Run_out.txt";
    string output_file = "LatticeCurrent.svg";
    double threshold = 1e-12;
    TotalMode total_mode = TOTAL_MODE_TOTAL;

    if (argc >= 4) {
        sites_per_unit_cell = atoi(argv[3]);
    }
    if (argc >= 5) {
        input_file = argv[4];
    }
    if (argc >= 6) {
        output_file = argv[5];
    }
    if (argc >= 7) {
        threshold = atof(argv[6]);
    }
    if (argc >= 8) {
        string mode_arg = argv[7];
        if (mode_arg == "total" || mode_arg == "Total" || mode_arg == "Total for set") {
            total_mode = TOTAL_MODE_TOTAL;
        }
        else if (mode_arg == "quantum" || mode_arg == "Quantum" || mode_arg == "Total quantum for set" || mode_arg == "Total_Quantum") {
            total_mode = TOTAL_MODE_QUANTUM;
        }
        else {
            cerr << "Invalid total_mode: " << mode_arg << endl;
            cerr << "Allowed values: total | quantum" << endl;
            return 1;
        }
    }

    if (Lx <= 0 || Ly <= 0 || sites_per_unit_cell <= 0) {
        cerr << "Lx, Ly, and sites_per_unit_cell must be positive integers." << endl;
        return 1;
    }

    vector<BondCurrent> bonds;
    int max_site = -1;
    if (!ParseRunOut(input_file, bonds, max_site, total_mode)) {
        return 1;
    }

    if (bonds.size() == 0) {
        if (total_mode == TOTAL_MODE_TOTAL) {
            cerr << "No valid (Set, sites, Total for set) entries found in " << input_file << endl;
        }
        else {
            cerr << "No valid (Set, sites, Total quantum for set) entries found in " << input_file << endl;
        }
        return 1;
    }

    int total_sites_from_output = max_site + 1;
    int expected_sites = Lx * Ly * sites_per_unit_cell;

    if (expected_sites != total_sites_from_output) {
        cerr << "Site-count mismatch: expected Lx*Ly*sites_per_unit_cell = " << expected_sites
             << ", but parsed max site index implies " << total_sites_from_output << " sites." << endl;
        cerr << "Parsed max site index = " << max_site << endl;
        return 1;
    }

    string total_label = (total_mode == TOTAL_MODE_TOTAL) ? "Total for set" : "Total quantum for set";

    if (!WriteSvg(output_file, bonds, expected_sites, Lx, Ly, sites_per_unit_cell, threshold, total_label)) {
        return 1;
    }

    cout << "Parsed sets: " << bonds.size() << endl;
    cout << "Output written: " << output_file << endl;
    cout << "Using: " << total_label << endl;
    cout << "Rule used: current>0 => m->n ; current<0 => n->m" << endl;

    return 0;
}
