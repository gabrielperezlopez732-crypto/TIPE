#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <random>
#include <chrono>
#include <iomanip>
#include "hashtable.h"
// #include "visualize_sprouts.h"


using namespace std;


#define M_PI 3.14159265358979323846
#define INFINITE 1000000

typedef std::string string;

typedef std::vector<int> vector_int;

typedef std::tuple<int, int> frontiere;

typedef std::tuple<int, int, int, vector<frontiere>> coup;




int int_of_string(string s) {
    return (int(s[0]) - int('0'))*10 + (int(s[1]) - int('0'));
}

string int_to_string(int x) {
    if (x >= 10)
        return std::to_string(x);
    else {
        string result = "0";
        result += std::to_string(x);
        return result;
    }
}

bool delimite_fr(string x) {
    return x == "++";
}

bool delimite_rg(string x) {
    return x == "**";
}

struct graphe {
    string graphe;
    vector_int regions; // indice de chaque région dans la chaîne
    vector<vector<frontiere>> frontieres_region; // Liste des frontieres dans chaque région
    vector<int> nb_frontieres; // Nombre de frontieres dans chaque région
    int nr; // nombre de régions
    int n; // nombre de points
    size_t taille;
    int* degres;
    bool** morts;
    tuple<frontiere, int>** frontieres;
    int*** accessibles;
    vector<int> vie_region;
};


typedef struct graphe graphe;

typedef tuple<int, string, bool> entree;


string extrait(graphe* g, int k) {
    if (k < 0 || k + 2 > (int)g->graphe.size()) {
            std::cerr << "extrait: invalid range k=" << k << " size=" << g->graphe.size() << std::endl;
            throw std::out_of_range("extrait: invalid range");
    }
    return g->graphe.substr(k, 2);
}

static string safe_substr(const string& src, size_t pos, size_t len) {
        if (pos > src.size()) {
            std::string snippet = src.substr(0, std::min((size_t)40, src.size()));
            std::cerr << "safe_substr: pos > src.size() pos=" << pos << " size=" << src.size() \
                      << " len=" << len << " src_snippet=\"" << snippet << "\"" << std::endl;
            throw std::out_of_range("safe_substr: pos > size");
        }
        size_t maxlen = len;
        if (pos + len > src.size() || len == string::npos) {
                maxlen = src.size() - pos;
        }
        return src.substr(pos, maxlen);
}

string connecter_string(string fi, string fj, int i, int j, int nouv) { 
    string f;
    
    f += safe_substr(fi, 0, i+2);

    f += int_to_string(nouv);


    f += safe_substr(fj, j, string::npos);     
    f += safe_substr(fj, 0, j);                

    f += int_to_string(nouv);

    f += safe_substr(fi, i+2, string::npos);
    
    return f;
}



std::tuple<frontiere, int> trouve_frontiere(graphe* g, int r, int i) {
    return g->frontieres[r][i];
}

void update_region_boundaries(graphe* g, int r) {
    if (r < 0 || r >= g->nr) return;

    int region_start = g->regions[r];
    int region_end = (r + 1 < (int)g->nr) ? g->regions[r + 1] : (int)g->taille;

    for (int k = region_start; k < region_end; k += 2) {
        if (k + 1 >= (int)g->taille) break;
        string x = extrait(g, k);
        if (x.length() == 2 && !delimite_fr(x) && !delimite_rg(x)) {
            int v = int_of_string(x);
            if (v >= 0 && v < g->n) {
                // Cherche le ++ (ou **) qui ouvre la face à gauche de k
                int fi1 = k - 2;
                while (fi1 >= 0 &&
                       !delimite_fr(extrait(g, fi1)) &&
                       !delimite_rg(extrait(g, fi1)))
                    fi1 -= 2;

                // Cherche le ++ (ou **) qui ferme la face à droite de k
                int fi2 = k + 2;
                while (fi2 < (int)g->taille &&
                       !delimite_fr(extrait(g, fi2)) &&
                       !delimite_rg(extrait(g, fi2)))
                    fi2 += 2;

                g->frontieres[r][v] = make_tuple(make_tuple(fi1, fi2), k);
            }
        }
    }
}

void rebuild_region_frontieres(graphe* g) {
    g->frontieres_region.clear();
    g->frontieres_region.resize(g->nr);
    g->nb_frontieres.assign(g->nr, 0);

    for (int r = 0; r < g->nr; ++r) {
        int region_start = g->regions[r];
        int region_end = (r + 1 < (int)g->nr) ? g->regions[r+1] : (int)g->taille;
        int k = region_start;
        while (k < region_end) {
            string tok = extrait(g, k);
            if (tok == "++") {
                int start = k;
                int kk = k + 2;
                while (kk < region_end) {
                    string tok2 = extrait(g, kk);
                    if (tok2 == "++" || tok2 == "**")
                        break;
                    kk += 2;
                }
                if (kk >= region_end)
                    break;
                g->frontieres_region[r].push_back(make_tuple(start, kk));
                g->nb_frontieres[r]++;
                k = kk;
            } else {
                k += 2;
            }
        }
    }
}

void rebuild_regions(graphe* g) {
    g->regions.clear();
    bool in_region = false;
    int N = (int)g->graphe.size();
    for (int k = 0; k + 1 < N; k += 2) {
        string tok = extrait(g, k);
        if (tok == "**") {
            if (!in_region) {
                g->regions.push_back(k);
                in_region = true;
            } else {
                in_region = false;
            }
        }
    }
    g->nr = (int)g->regions.size();
}

int compute_region_life_from_frontieres(graphe* g, int r) {
    vector<char> seen(g->n, 0);
    int life = 0;
    for (const auto& f : g->frontieres_region[r]) {
        int start = get<0>(f);
        int end = get<1>(f);
        for (int k = start; k < end; k += 2) {
            string tok = extrait(g, k);
            if (!delimite_fr(tok) && !delimite_rg(tok)) {
                int v = int_of_string(tok);
                if (v >= 0 && v < g->n && !seen[v]) {
                    seen[v] = 1;
                    life += max(0, 3 - g->degres[v]);
                }
            }
        }
    }
    return life;
}

string extraire (string s, int i, int j) {
    if (i < 0 || j < i || j > (int)s.size()) {
        std::cerr << "extraire: invalid indices i=" << i << " j=" << j << " size=" << s.size() << std::endl;
        throw std::out_of_range("extra ire: invalid indices");
    }
    return s.substr(i, j - i);
}

string couper_string(string s, int i, int j) {
    if (i < 0 || j < i || j > (int)s.size()) {
        std::cerr << "couper_string: invalid indices i=" << i << " j=" << j << " size=" << s.size() << std::endl;
        throw std::out_of_range("couper_string: invalid indices");
    }
    string s1 = s.substr(j);
    string s2 = s.substr(0, i);
    s2 += s1;
    return s2;
    
}

string couper_graphe(graphe* g, int i, int j) {
    // Renvoie la partie de g->graphe entre i (inclus) et j (exclus) 
    // Supprime cette partie du graphe
    if (i < 0 || j < i || j > (int)g->taille) {
        std::cerr << "couper_graphe: invalid indices i=" << i << " j=" << j << " size=" << g->taille << std::endl;
        throw std::out_of_range("couper_graphe: invalid indices");
    }
    string s = g->graphe.substr(i, j - i);
    string s1 = g->graphe.substr(0, i) + g->graphe.substr(j);
    g->graphe = s1;
    return s;
}



vector<int> trouver_regions(graphe* g, int i) { 
    // Renvoit toutes les régions 
    // auxquelles i appartient
    vector<int> regions;
    int r = 0;
    int N = g->taille;
    for (int ind = 0; ind < N - 1; ind+=2) {
        string x = extrait(g, ind);
        if (delimite_rg(x))
            r = ind;
        if (!delimite_fr(x)) {
            if (int_of_string(x) == i) {
                regions.push_back(r);
            }
        }
    }
    return regions;
}

int count(string s, int i) {
    int somme = 0;
    int N = s.size();
    for (int k = 0; k < N - 1; k+=2) {
        if (k < 0 || k + 2 > (int)s.size()) continue;
        string x = s.substr(k, 2);
        if (!delimite_fr(x) && !delimite_rg(x)) {
            if (int_of_string(x) == i)
                somme++;
        }
    }
    return somme;
}


static bool frontier_contains_vertex(const string& graph,
    const frontiere& f,
    int vertex,
    int& occurrence_pos)
{
    int start = get<0>(f);
    int end = get<1>(f);
    for (int k = start + 2; k < end; k += 2) {
        if (k + 1 >= (int)graph.size())
            break;
        string tok = graph.substr(k, 2);
        if (!delimite_fr(tok) && !delimite_rg(tok) && int_of_string(tok) == vertex) {
            occurrence_pos = k;
            return true;
        }
    }
    return false;
}

vector<int> voisins(graphe* g, int i) {
    vector<int> regions = trouver_regions(g, i);
    vector<int> v;
    int N = regions.size();
    for (int idx = 0; idx < N; idx++) {
        int r = regions[idx];
        int region_start = g->regions[r];
        int region_end = (r + 1 < g->nr) ? g->regions[r + 1] : g->taille;
        for (int k = region_start + 2; k < region_end - 2; k += 2) {
            string x = extrait(g, k);
            if (!delimite_fr(x) && !delimite_rg(x) && int_of_string(x) == i) {
                if (k - 2 >= region_start) {
                    string xm = extrait(g, k - 2);
                    if (!delimite_fr(xm) && !delimite_rg(xm)) {
                        v.push_back(int_of_string(xm));
                    }
                }
                if (k + 2 < region_end) {
                    string xp = extrait(g, k + 2);
                    if (!delimite_fr(xp) && !delimite_rg(xp)) {
                        v.push_back(int_of_string(xp));
                    }
                }
            }
        }
    }
    return v;
}

int degre(graphe* g, int s) {
    return (voisins(g, s)).size();
}

char get(graphe* g, int k) {
    return g->graphe[k] + g->graphe[k+1];
}

static void parse_planar_graph_string(const string& s,
    vector<vector<int>>& faces,
    vector<int>& face_regions)
{
    faces.clear();
    face_regions.clear();
    vector<int> current_face;
    int region = 0;
    int N = static_cast<int>(s.size());

    for (int k = 0; k + 1 < N; k += 2) {
        string tok = s.substr(k, 2);
        if (tok == "++") {
            if (!current_face.empty()) {
                faces.push_back(current_face);
                face_regions.push_back(region);
                current_face.clear();
            }
        } else if (tok == "**") {
            if (!current_face.empty()) {
                faces.push_back(current_face);
                face_regions.push_back(region);
                current_face.clear();
            }
            region++;
        } else {
            int v = stoi(tok);
            current_face.push_back(v);
        }
    }

    if (!current_face.empty()) {
        faces.push_back(current_face);
        face_regions.push_back(region);
    }
}

static void build_graph_edges(
    const vector<vector<int>>& faces,
    set<pair<int,int>>& edges,
    set<int>& vertices)
{
    edges.clear();
    vertices.clear();
    for (const auto& face : faces) {
        int m = static_cast<int>(face.size());
        if (m < 2) {
            continue;
        }
        for (int i = 0; i < m; ++i) {
            int a = face[i];
            int b = face[(i + 1) % m];
            if (a == b) {
                continue;
            }
            int u = min(a, b);
            int v = max(a, b);
            edges.insert(make_pair(u, v));
            vertices.insert(a);
            vertices.insert(b);
        }
    }
}

static void layout_vertices(
    const vector<int>& vertex_ids,
    const set<pair<int,int>>& edges,
    vector<double>& xs,
    vector<double>& ys)
{
    int n = static_cast<int>(vertex_ids.size());
    xs.assign(n, 0.0);
    ys.assign(n, 0.0);
    if (n == 0) {
        return;
    }

    for (int i = 0; i < n; ++i) {
        double angle = 2.0 * M_PI * i / max(1, n);
        xs[i] = cos(angle);
        ys[i] = sin(angle);
    }

    map<int,int> index_of;
    for (int i = 0; i < n; ++i) {
        index_of[vertex_ids[i]] = i;
    }

    int iterations = 120;
    double area = 1.0;
    double k = sqrt(area / max(1, n));

    vector<double> disp_x(n, 0.0);
    vector<double> disp_y(n, 0.0);

    for (int it = 0; it < iterations; ++it) {
        fill(disp_x.begin(), disp_x.end(), 0.0);
        fill(disp_y.begin(), disp_y.end(), 0.0);

        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double dx = xs[i] - xs[j];
                double dy = ys[i] - ys[j];
                double dist = sqrt(dx * dx + dy * dy) + 1e-6;
                double rep = (k * k) / dist;
                disp_x[i] += dx / dist * rep;
                disp_y[i] += dy / dist * rep;
                disp_x[j] -= dx / dist * rep;
                disp_y[j] -= dy / dist * rep;
            }
        }

        for (const auto& edge : edges) {
            int u = index_of.at(edge.first);
            int v = index_of.at(edge.second);
            double dx = xs[u] - xs[v];
            double dy = ys[u] - ys[v];
            double dist = sqrt(dx * dx + dy * dy) + 1e-6;
            double attr = (dist * dist) / k;
            disp_x[u] -= dx / dist * attr;
            disp_y[u] -= dy / dist * attr;
            disp_x[v] += dx / dist * attr;
            disp_y[v] += dy / dist * attr;
        }

        double temperature = 0.1 * (1.0 - static_cast<double>(it) / iterations);
        for (int i = 0; i < n; ++i) {
            double dist = sqrt(disp_x[i] * disp_x[i] + disp_y[i] * disp_y[i]);
            if (dist > 1e-6) {
                xs[i] += disp_x[i] / dist * min(dist, temperature);
                ys[i] += disp_y[i] / dist * min(dist, temperature);
            }
            xs[i] = min(max(xs[i], -2.0), 2.0);
            ys[i] = min(max(ys[i], -2.0), 2.0);
        }
    }
}

void afficher_graphe_planar(const string& s)
{
    vector<vector<int>> faces;
    vector<int> face_regions;
    parse_planar_graph_string(s, faces, face_regions);

    set<pair<int,int>> edges;
    set<int> vertices;
    build_graph_edges(faces, edges, vertices);

    vector<int> vertex_ids(vertices.begin(), vertices.end());
    vector<double> xs;
    vector<double> ys;
    layout_vertices(vertex_ids, edges, xs, ys);
    map<int,int> index_of;
    for (int i = 0; i < static_cast<int>(vertex_ids.size()); ++i) {
        index_of[vertex_ids[i]] = i;
    }

    cout << "Planar graph representation:" << endl;
    cout << "Vertices (" << vertex_ids.size() << "):" << endl;
    for (int i = 0; i < static_cast<int>(vertex_ids.size()); ++i) {
        cout << "  " << vertex_ids[i] << " -> (" << xs[i] << ", " << ys[i] << ")" << endl;
    }

    cout << "Edges (" << edges.size() << "):" << endl;
    for (const auto& edge : edges) {
        cout << "  " << edge.first << " - " << edge.second << endl;
    }

    cout << "Faces (" << faces.size() << "):" << endl;
    for (size_t i = 0; i < faces.size(); ++i) {
        cout << "  Face " << i << " (region " << face_regions[i] << "):";
        for (int vertex : faces[i]) {
            cout << " " << vertex;
        }
        cout << endl;
    }
}

void free_matrice(int** t, int l) {
    for (int i = 0; i < l; i++)
        free(t[i]);
    free(t);
}

int* degres_complet(graphe* g) {
    int n = g->n;
    int** vus = (int**)malloc(n*sizeof(int*));
    // vus[i][j] = true ssi on a compté l'arête (i, j) dans le sens
    // du parcours du graphe
    for (int i = 0 ; i < n; i++) 
        vus[i] = (int*)malloc(n*sizeof(int));
    int* deg = (int*)malloc(n*sizeof(int));
    int l = g->taille;
    for (int k = 2; k < l - 2; k+=2) {
        string x = extrait(g, k);
        string xm = extrait(g, k-2);
        string xp = extrait(g, k+2);
        if (!delimite_fr(x) && !delimite_rg(x)) {
            int i = int_of_string(x);
            if (!delimite_fr(xm) && !delimite_rg(xm)) {
                int im = int_of_string(xm);
                if (!vus[i][im]) {
                    vus[i][im] = true;
                    deg[i]+=1;
                    deg[im]+=1;
                }
            }
            if (!delimite_fr(xp) && !delimite_rg(xp)) {
                int ip = int_of_string(xp);
                if (!vus[i][ip]) {
                    vus[i][ip] = true;
                    deg[i]+=1;
                    deg[ip]+=1;
                }
            }
        }
    }
    free_matrice(vus, n);
    return(deg);
}


tuple<int, int, int, int> sommets_degres(graphe* g) {
    tuple<int, int, int, int> v = {0, 0, 0, 0};
    for (int i = 0; i < g->n; i++) {
        int d = degre(g, i);
        if (d == 0)
            get<0>(v) += 1;
        if (d == 1)
            get<1>(v) += 1;
        if (d == 2)
            get<2>(v) += 1;
        if (d == 3)
            get<3>(v) += 1;
    }
    return v;
}


vector<int> occurences_sommet (int i, string s) {
    vector<int> occ;
    int n = s.length();
    for (int j = 0; j < n - 2; j += 2) {
        if (j < 0 || j + 2 > (int)s.size()) continue;
        string x = s.substr(j, 2);
        if (!delimite_fr(x) && !delimite_rg(x)) {
            if (int_of_string(x) == i)
            occ.push_back(j);
        }
    }
    return occ;
}

int premiere_occ(int i, string s) {
    vector<int> occ;
    int n = s.length();
    for (int j = 0; j < n - 2; j += 2) {
        if (j < 0 || j + 2 > (int)s.size()) continue;
        string x = s.substr(j, 2);
        if (!delimite_fr(x) && !delimite_rg(x)) {
            if (int_of_string(x) == i)
                return j;
        }
    }
    return -1;
}

int derniere_occ(int i, string s) {
    int n = s.length();
    for (int j = n - 2; j >= 0; j -= 2) {
        if (j < 0 || j + 2 > (int)s.size()) continue;
        string x = s.substr(j, 2);
        if (!delimite_fr(x) && !delimite_rg(x)) {
            if (int_of_string(x) == i)
                return j;
        }
    }
    return -1;
}


void allocate_region_frontier(graphe* g, int new_region_index);
void update_region_boundaries(graphe* g, int r);
vector<int> liste_sommets_region(graphe* g, int r);

void connecter(graphe* g, int r, int i, int j, bool* region_modif, vector<frontiere> frontieres) {
    assert (g->degres[i] < 3 && "i est mort (degre >= 3)");
    assert (g->degres[j] < 3 && "j est mort (degre >= 3)");
    tuple<frontiere, int> ti = trouve_frontiere(g, r, i);
    tuple<frontiere, int> tj = trouve_frontiere(g, r, j);
    frontiere fi = get<0>(ti);
    frontiere fj = get<0>(tj);

    int indi = get<1>(ti);
    int indj = get<1>(tj);

    int fi1 = get<0>(fi);
    int fi2 = get<1>(fi);
    int fj1 = get<0>(fj);
    int fj2 = get<1>(fj);

    int old_life = g->vie_region[r];
    int old_nr = g->nr;
    bool new_region_created = false;

    auto find_selected_frontiere = [&](int vertex) {
        int occurrence_pos = -1;
        for (const auto& f : frontieres) {
            if (frontier_contains_vertex(g->graphe, f, vertex, occurrence_pos)) {
                return make_tuple(f, occurrence_pos);
            }
        }
        return trouve_frontiere(g, r, vertex);
    };

    if (!frontieres.empty()) {
        ti = find_selected_frontiere(i);
        tj = find_selected_frontiere(j);
        fi = get<0>(ti);
        fj = get<0>(tj);
        fi1 = get<0>(fi);
        fi2 = get<1>(fi);
        fj1 = get<0>(fj);
        fj2 = get<1>(fj);
        indi = get<1>(ti);
        indj = get<1>(tj);
    }

    if (fi != fj) {
        string si, sj;
        int local_i = indi - fi1;
        int local_j = indj - fj1;
        if (fi1 > fj1) {
            si = couper_graphe(g, fi1, fi2);
            sj = couper_graphe(g, fj1, fj2);
        } else {
            sj = couper_graphe(g, fj1, fj2);
            si = couper_graphe(g, fi1, fi2);
        }
        string s;
        s += "++";
        s += connecter_string(si, sj, local_i, local_j, g->n);
        s += "++";

        int insert_pos = min(fi1, fj1);
        g->graphe.insert(insert_pos, s);
    }

    else {
        *region_modif = true;
        new_region_created = true;
        int ind1 = min(indi, indj);
        int ind2 = max(indi, indj);

        string si = couper_graphe(g, fi1, fi2);
        int local_ind1 = ind1 - fi1;
        int local_ind2 = ind2 - fi1;

        string s;
        s += safe_substr(si, 0, local_ind1+2);        
        s += int_to_string(g->n);
        s += safe_substr(si, local_ind2, string::npos);


        string s0;
        s0 += safe_substr(si, 0, local_ind2+2);         
        s0 += int_to_string(g->n);
        s0 += safe_substr(si, local_ind1, string::npos);

        g->graphe.insert(fi1, s0);
        string appended = string("**++") + s + "++**";
        g->graphe.append(appended);
    }

    g->taille = g->graphe.size();
    rebuild_regions(g);
    if (g->nr > old_nr) {
        g->frontieres = (tuple<frontiere, int>**)realloc(g->frontieres, g->nr * sizeof(tuple<frontiere, int>*));
        g->frontieres[g->nr - 1] = nullptr; // sera alloué ci-dessous avec les autres
    }

    // 1. Enregistrer le nouveau sommet et étendre g->degres
    g->degres = (int*)realloc(g->degres, (g->n + 1) * sizeof(int));
    g->degres[g->n] = 2;  // tout nouveau sommet a exactement degré 2
    g->n++;               // g->n incrémenté AVANT update_region_boundaries

    // 2. Redimensionner tous les tableaux de frontières pour inclure le nouveau sommet
    for (int rr = 0; rr < g->nr; rr++) {
        tuple<frontiere, int>* nf = (tuple<frontiere, int>*)malloc(g->n * sizeof(tuple<frontiere, int>));
        if (!nf) {
            throw std::bad_alloc();
        }
        if (g->frontieres[rr]) {
            for (int vi = 0; vi < g->n - 1; vi++)
                nf[vi] = g->frontieres[rr][vi];
            nf[g->n - 1] = make_tuple(make_tuple(-1, -1), -1);
            free(g->frontieres[rr]);
        } else {
            for (int vi = 0; vi < g->n; vi++)
                nf[vi] = make_tuple(make_tuple(-1, -1), -1);
        }
        g->frontieres[rr] = nf;
    }

    // 3. Recalculer les frontières de toutes les régions
    for (int rr = 0; rr < g->nr; rr++) {
        update_region_boundaries(g, rr);
    }

    g->taille = g->graphe.size();
    rebuild_region_frontieres(g);

    if (g->accessibles != nullptr) {
        int new_n = g->n;

        for (int rr = 0; rr < old_nr; rr++) {
            g->accessibles[rr] = (int**)realloc(g->accessibles[rr], new_n * sizeof(int*));
            for (int a = 0; a < new_n - 1; a++) {
                g->accessibles[rr][a] = (int*)realloc(g->accessibles[rr][a], new_n * sizeof(int));
                g->accessibles[rr][a][new_n - 1] = 0;
            }
            g->accessibles[rr][new_n - 1] = (int*)calloc(new_n, sizeof(int));
        }

        if (g->nr > old_nr) {
            g->accessibles = (int***)realloc(g->accessibles, g->nr * sizeof(int**));
            g->accessibles[g->nr - 1] = (int**)malloc(new_n * sizeof(int*));
            for (int a = 0; a < new_n; a++)
                g->accessibles[g->nr - 1][a] = (int*)calloc(new_n, sizeof(int));
        }

        for (int rr = 0; rr < g->nr; rr++) {
            for (int a = 0; a < new_n; a++)
                for (int b = 0; b < new_n; b++)
                    g->accessibles[rr][a][b] = 0;

            for (const auto& f : g->frontieres_region[rr]) {
                int start = get<0>(f);
                int end   = get<1>(f);
                vector<int> verts;
                for (int k = start + 2; k < end; k += 2) {
                    if (k + 1 >= (int)g->taille) break;
                    string tok = g->graphe.substr(k, 2);
                    if (!delimite_fr(tok) && !delimite_rg(tok)) {
                        int v = int_of_string(tok);
                        if (v >= 0 && v < new_n)
                            verts.push_back(v);
                    }
                }
                for (int a : verts)
                    for (int b : verts)
                        g->accessibles[rr][a][b] = 1;
            }
        }
    }


    g->vie_region.resize(g->nr, 0);
    int delta_i = (i == j) ? 2 : 1;
    int delta_j = (i == j) ? 0 : 1;
    for (int rr = 0; rr < g->nr; rr++) {
        vector<int> vertices = liste_sommets_region(g, rr);
        int total_lives = 0;
        for (int v : vertices) {
            int deg = g->degres[v];
            if (v == i) deg += delta_i;
            if (v == j && i != j) deg += delta_j;
            total_lives += max(0, 3 - deg);
        }
        g->vie_region[rr] = max(0, total_lives - 1);
    }
}

bool est_dans_region(graphe* g, int r, int k) {
    if (r >= g->nr || r < 0) {
        return false;
    }
    int taille = (int)g->taille;
    if (r == g->nr - 1) {
        return g->regions[r] <= k && k < taille - 2;
    }
    else if (r < g->nr - 1) {
        return g->regions[r] <= k && k < g->regions[r+1];
    }
    return false;
}
    
vector<int> liste_sommets_region(graphe* g, int r) {
    int* vus = (int*)malloc(g->n*sizeof(int));
    for (int i = 0; i < g->n; i++) {
        vus[i] = 0;  // false
    }
    vector<int> liste;
    
    int region_start = g->regions[r];
    int region_end = (r + 1 < (int)g->nr) ? g->regions[r + 1] : (int)g->taille;
    
    for (int k = region_start; k < region_end; k += 2) {
        if (k + 1 < (int)g->taille) {
            string x = extrait(g, k);
            if (x.length() == 2 && !delimite_fr(x) && !delimite_rg(x)) {
                int i = int_of_string(x);
                if (i >= 0 && i < g->n && !vus[i]) {
                    vus[i] = 1;  // true
                    liste.push_back(i);
                }
            }
        }
    }
    free(vus);
    return liste;
}

vector<vector<frontiere>> sous_ensembles(const vector<frontiere>& f) {
    // Renvoie tous les sous_ensembles de f
    vector<vector<frontiere>> sous_ensembles;
    int n = f.size();
    for (int mask = 1; mask < (1 << n); mask++) {
        vector<frontiere> subset;
        for (int j = 0; j < n; j++) {
            if ((mask >> j) & 1) {
                subset.push_back(f[j]);
            }
        }
        sous_ensembles.push_back(subset);
    }
    return sous_ensembles;
}

static bool vertex_in_region(graphe* g, int r, int i) {
    if (r < 0 || r >= g->nr || i < 0 || i >= g->n) {
        return false;
    }
    int region_start = g->regions[r];
    int region_end = (r + 1 < g->nr) ? g->regions[r + 1] : (int)g->taille;
    for (int k = region_start; k < region_end; k += 2) {
        if (k + 1 >= (int)g->taille) {
            break;
        }
        string tok = extrait(g, k);
        if (!delimite_fr(tok) && !delimite_rg(tok) && int_of_string(tok) == i) {
            return true;
        }
    }
    return false;
}

vector<coup> coups(graphe* g) {
    vector<coup> coups_liste;
    for (int r = 0; r < g->nr; r++) {
        vector<int> vertices = liste_sommets_region(g, r);
        int m = vertices.size();
        const vector<frontiere>& fronts = g->frontieres_region[r];
        for (int a = 0; a < m; a++) {
            for (int b = a; b < m; b++) {
                int i = vertices[a];
                int j = vertices[b];
                bool legal_degree = (i == j) ? (g->degres[i] <= 1) : (g->degres[i] < 3 && g->degres[j] < 3);
                if (!legal_degree) {
                    continue;
                }
                int occ = -1;
                int fi_idx = -1, fj_idx = -1;
                for (int f = 0; f < (int)fronts.size(); f++) {
                    if (fi_idx == -1 && frontier_contains_vertex(g->graphe, fronts[f], i, occ))
                        fi_idx = f;
                    if (fj_idx == -1 && frontier_contains_vertex(g->graphe, fronts[f], j, occ))
                        fj_idx = f;
                }
                if (fi_idx == -1 || fj_idx == -1 || fi_idx == fj_idx) {

                    for (const auto& subset : sous_ensembles(fronts)) {
                        coups_liste.push_back(make_tuple(r, i, j, subset));
                    }
                } else {
                    coups_liste.push_back(make_tuple(r, i, j, vector<frontiere>{}));
                }
            }
        }
    }
    return coups_liste;
}


bool aucun_coup_region(graphe* g, int r) {
    return g->vie_region[r] == 0;
}

graphe* init_graphe(int n) {
    // Sommets numérotés de 0 à n-1
    graphe* g = new graphe();
    g->n = n;
    g->nr = 1;
    g->regions = {0};
    g->graphe = "****";
    for (int k = 0; k < n; k++) {
        g->graphe.insert(2, "++" + int_to_string(k) + "++");
    } 
    g->degres = (int*)malloc(n*sizeof(int));
    for (int k = 0; k < n; k++) {
        g->degres[k] = 0;
    }
    g->taille = 4+6*n;
    rebuild_region_frontieres(g);
    g->morts = (bool**)malloc(n*sizeof(bool*));
    for (int k = 0; k < n; k++) {
        g->morts[k] = (bool*)malloc(n*sizeof(bool));
        for (int l = 0; l < n; l++) {
            g->morts[k][l] = false;
        }
    }
    g->frontieres = (tuple<frontiere, int>**)malloc(g->nr*sizeof(tuple<frontiere, int>*));
    // g->frontiere[r][i] est la frontière de r à laquelle i appartient
    // et tuple<-1, -1> si i n'appartient à aucune frontière de r
    g->frontieres[0] = (tuple<frontiere, int>*)malloc(n*sizeof(tuple<frontiere, int>));
    for (int k = 0; k < n; k++) {
        g->frontieres[0][n-1-k] = make_tuple(make_tuple(2+6*k, 6+6*k), 4+6*k);
    }
    g->accessibles = (int***)malloc(g->nr * sizeof(int*));
    for (int r = 0; r < g->nr; r++) {
        g->accessibles[r] = (int**)malloc(n * sizeof(int*));
        for (int i = 0; i < n; i++) {
            g->accessibles[r][i] = (int*)malloc(n * sizeof(int));
            for (int j = 0; j < n; j++) {
                g->accessibles[r][i][j] = 1; // true: i et j sont accessibles dans la région r
            }
        }
    }
    g->vie_region = {3*n};
    return g;
}

void update_mort(graphe* g, int i, int r) {
    if (g->degres[i] >= 3) {
        g->morts[i][r] = true;
    }
}

void allocate_region_frontier(graphe* g, int new_region_index) {
    g->frontieres[new_region_index] = (tuple<frontiere, int>*)malloc(g->n*sizeof(tuple<frontiere, int>));
    for (int k = 0; k < g->n; k++) {
        g->frontieres[new_region_index][k] = make_tuple(make_tuple(-1, -1), -1);
    }
}


graphe maj(graphe g, coup c) {
    // On suppose que le coup (r, i, j) est valide pour g
    int r = get<0>(c);
    int i = get<1>(c);
    int j = get<2>(c);


    int* new_degres = (int*)malloc(g.n * sizeof(int));
    memcpy(new_degres, g.degres, g.n * sizeof(int));
    g.degres = new_degres;

    tuple<frontiere, int>** new_frontieres = (tuple<frontiere, int>**)malloc(g.nr * sizeof(tuple<frontiere, int>*));
    for (int rr = 0; rr < g.nr; rr++) {
        new_frontieres[rr] = (tuple<frontiere, int>*)malloc(g.n * sizeof(tuple<frontiere, int>));
        memcpy(new_frontieres[rr], g.frontieres[rr], g.n * sizeof(tuple<frontiere, int>));
    }
    g.frontieres = new_frontieres;
    g.morts = nullptr;
    g.accessibles = nullptr;

    bool b = false;
    connecter(&g, r, i, j, &b, get<3>(c));
    if (i == j) {
        g.degres[i] += 2;
    }
    else {
        g.degres[i]++;
        g.degres[j]++;
    }
    if (b) {

        int new_r = g.nr - 1;
        if (!g.frontieres_region[new_r].empty()) {
            g.frontieres_region[new_r].erase(g.frontieres_region[new_r].begin());
            g.nb_frontieres[new_r]--;
        }
    }
    return g;
}

void maj_statique(graphe* g, coup c) {
    // On suppose que le coup (r, i, j) est valide pour g
    int r = get<0>(c);
    int i = get<1>(c);
    int j = get<2>(c);
    vector<frontiere> frontieres = get<3>(c);
    bool b = false;
    connecter(g, r, i, j, &b, frontieres);
    if (i == j) {
        g->degres[i] += 2;
    }
    else {
        g->degres[i]++;
        g->degres[j]++;
    }
    if (b) {
        int new_r = g->nr - 1;
        if (!g->frontieres_region[new_r].empty()) {
            g->frontieres_region[new_r].erase(g->frontieres_region[new_r].begin());
            g->nb_frontieres[new_r]--;
        }
    }
}


int minimum(int x, int y) {
    if (x < y)
        return x;
    else
        return y;
}

int maximum(int x, int y) {
    if (x > y)
        return x;
    else
        return y;
}



float proportion_coups(graphe* g, int l) {
    int total_dead_regions = 0;
    for (int r = 0; r < g->nr; r++) {
        if (aucun_coup_region(g, r)) {
            total_dead_regions += 1;
        }
    }
    if (l <= 0) {
        return 0.0f;
    }
    return (float)total_dead_regions / (float)l;
}

float g_scaling_factor = 10.0f;

float proportion(graphe* g, int l) {
    float ratio = proportion_coups(g, l);
    float value = fabs((1.0f - float(l) / g_scaling_factor) * ratio);
    float result = std::min(1.0f, value);

    return std::max(0.1f, result);
}



vector<int> arr_to_vect (int* t, int len) {
    vector<int> v;
    for (int i = 0; i < len; i++) {
        v.push_back(t[i]);
    }
    return v;
}

bool coup_legal(graphe* g, coup c) {
    int r = get<0>(c);
    int i = get<1>(c);
    int j = get<2>(c);
    vector<frontiere> frontieres = get<3>(c);
    if (r < 0 || r >= g->nr || i < 0 || i >= g->n || j < 0 || j >= g->n) {
        return false;
    }
    bool legal_degree = (i == j) ? (g->degres[i] <= 1) : (g->degres[i] < 3 && g->degres[j] < 3);
    if (!legal_degree) {
        return false;
    }
    if (!vertex_in_region(g, r, i) || !vertex_in_region(g, r, j)) {
        return false;
    }
    for (const auto& f : frontieres) {
        int start = get<0>(f);
        int end = get<1>(f);
        if (!est_dans_region(g, r, start) || !est_dans_region(g, r, end)) {
            return false;
        }
    }
    return true;
}

int compte_sommets_vivants(graphe* g) {
    int count = 0;
    for (int i = 0; i < g->n; i++) {
        if (g->degres[i] < 3) {
            count++;
        }
    }
    return count;
}

int score_naif(graphe* g, int joueur) {
    return 0;
}

int score_regions_max(graphe* g, int joueur) {
    int score = 0;
    for (int r = 0 ; r < g->nr; r++) {
        if (aucun_coup_region(g, r)) {}
        score += 1;
    }
    return score;
}

int score_regions_min(graphe* g, int joueur) {
    int score = 0;
    for (int r = 0 ; r < g->nr; r++) {
        if (aucun_coup_region(g, r)) {}
        score -= 1;
    }
    return score;
}

int score_regions_vivants(graphe* g, int joueur) {
    int score = 0;
    for (int r = 0 ; r < g->nr; r++) {
        if (aucun_coup_region(g, r)) {
            if (joueur == 1)
                score -= 5;   
            else
                score += 5;  
        }
    }
    score += compte_sommets_vivants(g) * (joueur == 1 ? 1 : -1);
    return score;
}

int score_sommets_parite(graphe * g, int joueur) {
    int score = 0;
    for (int r = 0 ; r < g->nr; r++) {
        if (g->vie_region[r] > 0) {
            int nb_vivants = g->vie_region[r];
            if (nb_vivants % 2 == 1) {
                score += (joueur == 1) ? 3 : -3;
            } else {
                score += (joueur == 1) ? -3 : 3;
            }
        }
    }
    return score;
}




int eval(graphe* g, vector<coup> coups, int joueur, int g_eval_type) {
    int l = coups.size();
    if (l == 0) {
        if (joueur == 1)
            return -INFINITE; 
        else
            return INFINITE;   
    }
    int score;
    switch(g_eval_type) {
        case 0: score = score_naif(g, joueur); break;
        case 1: score = score_regions_max(g, joueur); break;
        case 2: score = score_regions_min(g, joueur); break;
        case 3: score = score_sommets_parite(g, joueur); break;
        default: score = score_naif(g, joueur);
    }
    return score;
}




int calcule_profondeur(graphe* g, vector<coup> coups_liste, int max_profondeur, int c) {
    int l = coups_liste.size();
    if (l == 0)
        return 0;
    return min(max_profondeur, 1 + c * (int)log2(l));
}

int alpha_beta(graphe* g, int alpha, int beta, int joueur, int profondeur, HashTable* h, int g_eval_type) {
    // L'IA est max = 1, le joueur est min = 0
    if (h->contains(g->graphe)) {
        return h->get(g->graphe);
    }
    vector<coup> coups_liste = coups(g);
    const int MAX_BRANCH = 10;
    if ((int)coups_liste.size() > MAX_BRANCH)
        coups_liste.resize(MAX_BRANCH);
    int l = coups_liste.size();
    int v;
    if (profondeur == 0 || l == 0) {
        v = eval(g, coups_liste, joueur, g_eval_type);
        h->insert(g->graphe, v);
        return v;
    }

    if (joueur == 1) {
        v = -INFINITE;
        for (int k = 0 ; k < l ; k++) {
            graphe g_nouv = maj(*g, coups_liste[k]);
            v = maximum(v, alpha_beta(&g_nouv, alpha, beta, 1-joueur, profondeur-1, h   , g_eval_type));
            if (v >= beta) {
                h->insert(g->graphe, v);
                return v;
            }
            alpha = maximum(alpha, v);
        }
    }
    else {
        v = INFINITE;
        for (int k = 0 ; k < l ; k++) {
            graphe g_nouv = maj(*g, coups_liste[k]);
            v = minimum(v, alpha_beta(&g_nouv, alpha, beta, 1-joueur, profondeur-1, h, g_eval_type));
            if (v <= alpha) {
                h->insert(g->graphe, v);
                return v;
            }
            beta = minimum(beta, v);
        }
    }
    h->insert(g->graphe, v);
    return v;
}


coup jeu_naif(vector<coup> coups_liste, int l) {
    int n = rand() % l;
    return coups_liste[n];
}


int main(void) {
    int g_eval_type = 3;  // score_sommets_parite
    int max_profondeur = 4;
    int c = 5;

    cout << "\n========== score_sommets_parite AI vs AI (n=5, 5 games) ==========" << endl;
    cout << "====================================\n" << endl;
    cout << "n | P0 | P1 | RUNTIME (s)" << endl;
    cout << "--|----|----|------------" << endl;

    for (int n = 5; n <= 5; n++) {
        auto debut = std::chrono::high_resolution_clock::now();
        int p0_wins = 0, p1_wins = 0;

        for (int game = 0; game < 5; game++) {
            HashTable* h = new HashTable(1000);
            graphe* g = init_graphe(n);
            std::mt19937 rng(game * 13 + 7);
            int player = 1;
            int moves = 0;

            while (moves < 200) {
                vector<coup> moves_list = coups(g);
                if (moves_list.empty()) {
                    if (player == 0) p0_wins++;
                    else p1_wins++;
                    break;
                }

                vector<coup> search_moves = moves_list;
                int depth = calcule_profondeur(g, search_moves, max_profondeur, c);
                float fraction = proportion(g, (int)moves_list.size());
                int sample_size = std::max(1, (int)std::ceil(fraction * (int)moves_list.size()));
                std::shuffle(search_moves.begin(), search_moves.end(), rng);
                if (sample_size < (int)search_moves.size()) search_moves.resize(sample_size);

                int best_score = -INFINITE;
                coup chosen = search_moves[0];
                int opponent = 1 - player;

                for (const auto& m : search_moves) {
                    graphe g_new = maj(*g, m);
                    int score = alpha_beta(&g_new, -INFINITE, INFINITE, opponent, depth, h, g_eval_type);
                    if (score > best_score) {
                        best_score = score;
                        chosen = m;
                    }
                }

                maj_statique(g, chosen);
                player = 1 - player;
                moves++;
            }
            delete h;
        }

        auto fin = std::chrono::high_resolution_clock::now();
        double runtime = std::chrono::duration<double>(fin - debut).count();

        cout << n << " | " << std::setw(2) << p0_wins << " | " << std::setw(2) << p1_wins << " | "
             << std::fixed << std::setprecision(2) << runtime << endl;
    }

    cout << "\n========== DONE ==========" << endl;
    return 0;
}


