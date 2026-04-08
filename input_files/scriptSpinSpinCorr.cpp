#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace std;

struct Point {
    double x;
    double y;
};

struct Bond {
    int i;
    int j;
    double corr;
    double dist;
};

static string ToHexColor(int r, int g, int b) {
    const char *hex = "0123456789abcdef";
    string out = "#000000";
    out[1] = hex[(r >> 4) & 0xF];
    out[2] = hex[r & 0xF];
    out[3] = hex[(g >> 4) & 0xF];
    out[4] = hex[g & 0xF];
    out[5] = hex[(b >> 4) & 0xF];
    out[6] = hex[b & 0xF];
    return out;
}

static string CorrToColor(double corr) {
    const double cmax = 0.75;
    if (corr >= 0.0) {
        double t = corr / cmax;
        if (t > 1.0) { t = 1.0; }
        int r = static_cast<int>(255.0 * (1.0 - t) + 0.5);
        int g = static_cast<int>(255.0 * (1.0 - t) + 0.5);
        int b = 255;
        return ToHexColor(r, g, b);
    }
    else {
        double t = (-corr) / cmax;
        if (t > 1.0) { t = 1.0; }
        int r = 255;
        int g = static_cast<int>(255.0 * (1.0 - t) + 0.5);
        int b = static_cast<int>(255.0 * (1.0 - t) + 0.5);
        return ToHexColor(r, g, b);
    }
}

static bool ReadCorrelationMatrix(const string &filepath,
                                  const string &block_header,
                                  vector<vector<double> > &mat) {
    ifstream infile(filepath.c_str());
    if (!infile.is_open()) {
        cerr << "Unable to open input file: " << filepath << endl;
        return false;
    }

    string start_tag = "--------------<" + block_header + ">-------------------";
    string end_tag = "-------------------------------------------------------";

    regex pair_re("\\(([-+0-9.eE]+),([-+0-9.eE]+)\\)");
    string line;
    bool in_block = false;
    mat.clear();

    while (getline(infile, line)) {
        if (!in_block) {
            if (line.find(start_tag) != string::npos) {
                in_block = true;
            }
            continue;
        }

        if (line.find(end_tag) != string::npos) {
            break;
        }

        vector<double> row;
        for (sregex_iterator it(line.begin(), line.end(), pair_re), it_end; it != it_end; ++it) {
            double real_part = atof((*it)[1].str().c_str());
            row.push_back(real_part);
        }

        if (!row.empty()) {
            mat.push_back(row);
        }
    }

    if (mat.empty()) {
        cerr << "Could not find <" << block_header << "> block in file." << endl;
        return false;
    }

    size_t nrows = mat.size();
    for (size_t r = 0; r < nrows; r++) {
        if (mat[r].size() != nrows) {
              cerr << "Parsed matrix for <" << block_header << "> is not square. Rows=" << nrows
                 << " row[" << r << "] has " << mat[r].size() << " entries." << endl;
            return false;
        }
    }

    return true;
}

static bool BuildTotalSpinSpinCorr(const vector<vector<double> > &szsz,
                                   const vector<vector<double> > &splus_sminus,
                                   const vector<vector<double> > &sminus_splus,
                                   vector<vector<double> > &total_corr) {
    if (szsz.size() != splus_sminus.size() || szsz.size() != sminus_splus.size()) {
        cerr << "Matrix size mismatch among SzSz, SplusSminus, and SminusSplus blocks." << endl;
        return false;
    }

    int n = static_cast<int>(szsz.size());
    total_corr.assign(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        if (static_cast<int>(szsz[i].size()) != n ||
            static_cast<int>(splus_sminus[i].size()) != n ||
            static_cast<int>(sminus_splus[i].size()) != n) {
            cerr << "Matrix row-size mismatch while building total spin-spin correlation." << endl;
            return false;
        }

        for (int j = 0; j < n; j++) {
            total_corr[i][j] = szsz[i][j] + 0.5 * (sminus_splus[i][j] + splus_sminus[i][j]);
        }
    }

    return true;
}

static Point SiteToModelPoint_2site_uc(int site, int Lx, int Ly) {
    int a = site % 2;
    int cell_linear = site / 2;
    int ux = cell_linear % Lx;
    int uy = cell_linear / Lx;

    Point p;
    p.x = ux + ((a == 0) ? 0.5 : 0.0);
    p.y = uy + ((a == 0) ? 0.0 : 0.5);
    return p;
}

static Point PeriodicDelta(const Point &from, const Point &to, double width, double height) {
    Point d;
    d.x = to.x - from.x;
    d.y = to.y - from.y;

    if (d.x > 0.5 * width) { d.x -= width; }
    if (d.x < -0.5 * width) { d.x += width; }
    if (d.y > 0.5 * height) { d.y -= height; }
    if (d.y < -0.5 * height) { d.y += height; }

    return d;
}

static bool IsInsideBox(const Point &p, double xmin, double xmax, double ymin, double ymax) {
    return (p.x >= xmin && p.x <= xmax && p.y >= ymin && p.y <= ymax);
}

static Point WrapPointToBox(Point p, double xmin, double xmax, double ymin, double ymax, double width, double height) {
    while (p.x < xmin) { p.x += width; }
    while (p.x > xmax) { p.x -= width; }
    while (p.y < ymin) { p.y += height; }
    while (p.y > ymax) { p.y -= height; }
    return p;
}

static vector<Bond> BuildNearestNeighborBonds(const vector<vector<double> > &mat, int Lx, int Ly) {
    int nsites = static_cast<int>(mat.size());

    vector<Point> model_pos(nsites);
    for (int s = 0; s < nsites; s++) {
        model_pos[s] = SiteToModelPoint_2site_uc(s, Lx, Ly);
    }

    double width = static_cast<double>(Lx);
    double height = static_cast<double>(Ly);

    double min_dist = 1e18;
    for (int i = 0; i < nsites; i++) {
        for (int j = i + 1; j < nsites; j++) {
            Point d = PeriodicDelta(model_pos[i], model_pos[j], width, height);
            double dist = sqrt(d.x * d.x + d.y * d.y);
            if (dist > 1e-12 && dist < min_dist) {
                min_dist = dist;
            }
        }
    }

    double tol = 1e-7;
    vector<Bond> bonds;
    for (int i = 0; i < nsites; i++) {
        for (int j = i + 1; j < nsites; j++) {
            Point d = PeriodicDelta(model_pos[i], model_pos[j], width, height);
            double dist = sqrt(d.x * d.x + d.y * d.y);
            if (fabs(dist - min_dist) <= tol) {
                Bond b;
                b.i = i;
                b.j = j;
                b.corr = 0.5 * (mat[i][j] + mat[j][i]);
                b.dist = dist;
                bonds.push_back(b);
            }
        }
    }

    return bonds;
}

static vector<Bond> BuildNextNearestNeighborBonds(const vector<vector<double> > &mat, int Lx, int Ly) {
    int nsites = static_cast<int>(mat.size());

    vector<Point> model_pos(nsites);
    for (int s = 0; s < nsites; s++) {
        model_pos[s] = SiteToModelPoint_2site_uc(s, Lx, Ly);
    }

    double width = static_cast<double>(Lx);
    double height = static_cast<double>(Ly);

    // Collect all unique pairwise distances
    vector<double> all_dists;
    for (int i = 0; i < nsites; i++) {
        for (int j = i + 1; j < nsites; j++) {
            Point d = PeriodicDelta(model_pos[i], model_pos[j], width, height);
            double dist = sqrt(d.x * d.x + d.y * d.y);
            if (dist > 1e-12) {
                all_dists.push_back(dist);
            }
        }
    }

    if (all_dists.empty()) { return vector<Bond>(); }
    sort(all_dists.begin(), all_dists.end());

    // Find the shortest (NN) distance
    double min_dist = all_dists[0];
    // Find the second-shortest (NNN) distance
    double tol = 1e-7;
    double nnn_dist = -1.0;
    for (size_t k = 0; k < all_dists.size(); k++) {
        if (all_dists[k] > min_dist + tol) {
            nnn_dist = all_dists[k];
            break;
        }
    }
    if (nnn_dist < 0.0) { return vector<Bond>(); }

    vector<Bond> bonds;
    for (int i = 0; i < nsites; i++) {
        for (int j = i + 1; j < nsites; j++) {
            Point d = PeriodicDelta(model_pos[i], model_pos[j], width, height);
            double dist = sqrt(d.x * d.x + d.y * d.y);
            if (fabs(dist - nnn_dist) <= tol) {
                Bond b;
                b.i = i;
                b.j = j;
                b.corr = 0.5 * (mat[i][j] + mat[j][i]);
                b.dist = dist;
                bonds.push_back(b);
            }
        }
    }

    return bonds;
}

static bool WriteSvg(const string &outfile,
                     const vector<Bond> &nn_bonds,
                     const vector<Bond> &nnn_bonds,
                     int Lx,
                     int Ly,
                     int nsites,
                     double threshold) {

    ofstream svg(outfile.c_str());
    if (!svg.is_open()) {
        cerr << "Unable to open output file: " << outfile << endl;
        return false;
    }

    const double unit = 140.0;
    const double margin = 70.0;
    const double lattice_w = Lx * unit;
    const double lattice_h = Ly * unit;
    const double legend_panel_w = 260.0;
    const double canvas_w = lattice_w + 2.0 * margin + legend_panel_w;
    const double canvas_h = lattice_h + 2.0 * margin;

    vector<Point> site_model(nsites), site_screen(nsites);
    for (int s = 0; s < nsites; s++) {
        site_model[s] = SiteToModelPoint_2site_uc(s, Lx, Ly);
        site_screen[s].x = margin + site_model[s].x * unit;
        site_screen[s].y = margin + (Ly - site_model[s].y) * unit;
    }

    double max_abs_corr = 0.0;
    for (size_t b = 0; b < nn_bonds.size(); b++) {
        max_abs_corr = max(max_abs_corr, fabs(nn_bonds[b].corr));
    }
    for (size_t b = 0; b < nnn_bonds.size(); b++) {
        max_abs_corr = max(max_abs_corr, fabs(nnn_bonds[b].corr));
    }
    if (max_abs_corr == 0.0) {
        max_abs_corr = 1.0;
    }

    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << canvas_w
        << "\" height=\"" << canvas_h << "\" viewBox=\"0 0 " << canvas_w << " " << canvas_h << "\">\n";

    svg << "<rect x=\"0\" y=\"0\" width=\"" << canvas_w << "\" height=\"" << canvas_h << "\" fill=\"white\"/>\n";
    svg << "<text x=\"" << 0.5 * (lattice_w + 2.0 * margin) << "\" y=\"34\" text-anchor=\"middle\" font-size=\"22\" font-family=\"Arial\">NN &amp; NNN total spin-spin: SzSz + 0.5(S-S+ + S+S-)</text>\n";

    svg << "<rect x=\"" << margin << "\" y=\"" << margin << "\" width=\"" << lattice_w << "\" height=\"" << lattice_h
        << "\" fill=\"#fafafa\" stroke=\"#d1d5db\" stroke-width=\"1.2\"/>\n";

    double width_model = static_cast<double>(Lx);
    double height_model = static_cast<double>(Ly);

    double xmin = margin;
    double xmax = margin + lattice_w;
    double ymin = margin;
    double ymax = margin + lattice_h;

    // Lambda-like helper to draw a list of bonds with a given dash pattern
    struct DrawBonds {
        static void draw(ofstream &svg, const vector<Bond> &bonds,
                         const vector<Point> &site_screen, const vector<Point> &site_model,
                         double width_model, double height_model, double unit,
                         double xmin, double xmax, double ymin, double ymax,
                         double lattice_w, double lattice_h,
                         double max_abs_corr, double threshold,
                         const string &dash) {
            for (size_t b = 0; b < bonds.size(); b++) {
                const Bond &bond = bonds[b];
                if (fabs(bond.corr) < threshold) { continue; }

                Point pi = site_screen[bond.i];
                Point d_model = PeriodicDelta(site_model[bond.i], site_model[bond.j], width_model, height_model);
                Point d_screen;
                d_screen.x = d_model.x * unit;
                d_screen.y = -d_model.y * unit;

                Point pj;
                pj.x = pi.x + d_screen.x;
                pj.y = pi.y + d_screen.y;

                double scaled = fabs(bond.corr) / max_abs_corr;
                double stroke_w = 1.0 + 7.0 * scaled;
                string color = CorrToColor(bond.corr);
                string dash_attr = dash.empty() ? "" : (" stroke-dasharray=\"" + dash + "\"");

                if (IsInsideBox(pj, xmin, xmax, ymin, ymax)) {
                    svg << "<line x1=\"" << pi.x << "\" y1=\"" << pi.y
                        << "\" x2=\"" << pj.x << "\" y2=\"" << pj.y
                        << "\" stroke=\"" << color << "\" stroke-width=\"" << stroke_w
                        << "\"" << dash_attr << ">\n";
                    svg << "<title>i=" << bond.i << " j=" << bond.j << " corr="
                        << scientific << setprecision(6) << bond.corr << "</title>\n";
                    svg << "</line>\n";
                } else {
                    double dx = pj.x - pi.x;
                    double dy = pj.y - pi.y;
                    double tx = 1e18, ty = 1e18;
                    if (dx > 0.0 && pj.x > xmax) { tx = (xmax - pi.x) / dx; }
                    if (dx < 0.0 && pj.x < xmin) { tx = (xmin - pi.x) / dx; }
                    if (dy > 0.0 && pj.y > ymax) { ty = (ymax - pi.y) / dy; }
                    if (dy < 0.0 && pj.y < ymin) { ty = (ymin - pi.y) / dy; }
                    double t = min(tx, ty);
                    if (!(t > 0.0 && t < 1.0)) { t = 0.5; }
                    Point p_cross;
                    p_cross.x = pi.x + t * dx;
                    p_cross.y = pi.y + t * dy;
                    Point p_cross_wrapped = WrapPointToBox(p_cross, xmin, xmax, ymin, ymax, lattice_w, lattice_h);
                    Point pj_wrapped = WrapPointToBox(pj, xmin, xmax, ymin, ymax, lattice_w, lattice_h);
                    svg << "<line x1=\"" << pi.x << "\" y1=\"" << pi.y
                        << "\" x2=\"" << p_cross.x << "\" y2=\"" << p_cross.y
                        << "\" stroke=\"" << color << "\" stroke-width=\"" << stroke_w
                        << "\"" << dash_attr << ">\n";
                    svg << "<title>i=" << bond.i << " j=" << bond.j << " corr="
                        << scientific << setprecision(6) << bond.corr << "</title>\n";
                    svg << "</line>\n";
                    svg << "<line x1=\"" << p_cross_wrapped.x << "\" y1=\"" << p_cross_wrapped.y
                        << "\" x2=\"" << pj_wrapped.x << "\" y2=\"" << pj_wrapped.y
                        << "\" stroke=\"" << color << "\" stroke-width=\"" << stroke_w
                        << "\"" << dash_attr << "/>\n";
                }
            }
        }
    };

    // Draw NNN first (behind NN)
    DrawBonds::draw(svg, nnn_bonds, site_screen, site_model,
                    width_model, height_model, unit,
                    xmin, xmax, ymin, ymax, lattice_w, lattice_h,
                    max_abs_corr, threshold, "");
    // Draw NN on top
    DrawBonds::draw(svg, nn_bonds, site_screen, site_model,
                    width_model, height_model, unit,
                    xmin, xmax, ymin, ymax, lattice_w, lattice_h,
                    max_abs_corr, threshold, "");

    for (int s = 0; s < nsites; s++) {
        svg << "<circle cx=\"" << site_screen[s].x << "\" cy=\"" << site_screen[s].y
            << "\" r=\"7.5\" fill=\"#111827\"/>\n";
        svg << "<text x=\"" << (site_screen[s].x + 10.0) << "\" y=\"" << (site_screen[s].y - 9.0)
            << "\" font-size=\"12\" font-family=\"Arial\" fill=\"#111827\">" << s << "</text>\n";
    }

    double lx = margin + lattice_w + 20.0;
    double ly = margin + 10.0;
    svg << "<rect x=\"" << lx << "\" y=\"" << ly << "\" width=\"220\" height=\"112\" fill=\"#ffffff\" stroke=\"#d1d5db\"/>\n";
    svg << "<text x=\"" << (lx + 14) << "\" y=\"" << (ly + 20) << "\" font-size=\"13\" font-family=\"Arial\">NN &amp; NNN bonds</text>\n";
    // Color scale reference lines
    svg << "<line x1=\"" << (lx + 16) << "\" y1=\"" << (ly + 34) << "\" x2=\"" << (lx + 84) << "\" y2=\"" << (ly + 34)
        << "\" stroke=\"#ff0000\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"" << (lx + 94) << "\" y=\"" << (ly + 39) << "\" font-size=\"12\" font-family=\"Arial\">corr = -0.75 (red)</text>\n";
    svg << "<line x1=\"" << (lx + 16) << "\" y1=\"" << (ly + 56) << "\" x2=\"" << (lx + 84) << "\" y2=\"" << (ly + 56)
        << "\" stroke=\"#ffffff\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"" << (lx + 94) << "\" y=\"" << (ly + 61) << "\" font-size=\"12\" font-family=\"Arial\">corr = 0 (white)</text>\n";
    svg << "<line x1=\"" << (lx + 16) << "\" y1=\"" << (ly + 78) << "\" x2=\"" << (lx + 84) << "\" y2=\"" << (ly + 78)
        << "\" stroke=\"#0000ff\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"" << (lx + 94) << "\" y=\"" << (ly + 83) << "\" font-size=\"12\" font-family=\"Arial\">corr = +0.75 (blue)</text>\n";


    svg << "</svg>\n";
    svg.close();
    return true;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        cout << "Usage: ./scriptSpinSpinCorr Lx Ly [Run_out.txt] [output.svg] [threshold]\n";
        cout << "Example: ./scriptSpinSpinCorr 3 2 Run_out.txt SpinSpinCorr.svg 1e-12\n";
        return 1;
    }

    int Lx = atoi(argv[1]);
    int Ly = atoi(argv[2]);
    string input_file = "Run_out.txt";
    string output_file = "SpinSpinCorr.svg";
    double threshold = 1e-12;

    if (argc >= 4) {
        input_file = argv[3];
    }
    if (argc >= 5) {
        output_file = argv[4];
    }
    if (argc >= 6) {
        threshold = atof(argv[5]);
    }

    if (Lx <= 0 || Ly <= 0) {
        cerr << "Lx and Ly must be positive integers." << endl;
        return 1;
    }

    vector<vector<double> > szsz, splus_sminus, sminus_splus, total_corr;
    if (!ReadCorrelationMatrix(input_file, "Sz[i].Sz[j]", szsz)) {
        return 1;
    }
    if (!ReadCorrelationMatrix(input_file, "Splus[i].Sminus[j]", splus_sminus)) {
        return 1;
    }
    if (!ReadCorrelationMatrix(input_file, "Sminus[i].Splus[j]", sminus_splus)) {
        return 1;
    }
    if (!BuildTotalSpinSpinCorr(szsz, splus_sminus, sminus_splus, total_corr)) {
        return 1;
    }

    int nsites = static_cast<int>(total_corr.size());
    int expected_sites = 2 * Lx * Ly;
    if (nsites != expected_sites) {
        cerr << "Site-count mismatch for 2-site unit-cell lattice: matrix has " << nsites
             << " but expected 2*Lx*Ly = " << expected_sites << endl;
        return 1;
    }

    vector<Bond> nn_bonds = BuildNearestNeighborBonds(total_corr, Lx, Ly);
    if (nn_bonds.empty()) {
        cerr << "No nearest-neighbor bonds identified." << endl;
        return 1;
    }

    vector<Bond> nnn_bonds = BuildNextNearestNeighborBonds(total_corr, Lx, Ly);

    if (!WriteSvg(output_file, nn_bonds, nnn_bonds, Lx, Ly, nsites, threshold)) {
        return 1;
    }

    cout << "Parsed matrices: SzSz, SplusSminus, SminusSplus (" << nsites << " x " << nsites << ")" << endl;
    cout << "Nearest-neighbor bonds plotted: " << nn_bonds.size() << endl;
    cout << "Next-nearest-neighbor bonds plotted: " << nnn_bonds.size() << endl;
    cout << "Using total corr = SzSz + 0.5*(SminusSplus + SplusSminus)" << endl;
    cout << "Output written: " << output_file << endl;

    return 0;
}
