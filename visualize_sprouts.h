#ifndef VISUALIZE_SPROUTS_H
#define VISUALIZE_SPROUTS_H

#include <map>
#include <set>
#include <string>
#include <vector>

// ── Types publics ─────────────────────────────────────────────────────────────

struct Vec2 { double x, y; };

struct Face {
    int region;
    int face_idx;
    std::vector<int> verts;   // sommets dans l'ordre RHS de la face
};

struct SproutsGraph {
    std::vector<int>                     vertices;
    std::set<std::pair<int,int>>         edges;
    std::set<int>                        self_loops;
    std::vector<Face>                    faces;
};

// Compute planar faces from the geometric layout
std::vector<Face> compute_planar_faces(const SproutsGraph& g, const std::map<int, Vec2>& pos);

// ── API publique ──────────────────────────────────────────────────────────────

// Parse la chaîne Sprouts et retourne le graphe extrait.
SproutsGraph parse_sprouts(const std::string& s);

// Calcule les positions 2-D des sommets (Fruchterman-Reingold).
// W, H : dimensions du canvas en pixels.
std::map<int, Vec2> compute_layout(const SproutsGraph& g, double W, double H);

// Génère un fichier SVG représentant le graphe.
// s        : chaîne Sprouts (ex. "**++02++++01++++00++**")
// filename : chemin du fichier de sortie (ex. "graph.svg")
void generate_svg(const std::string& s, const std::string& filename);

#endif // VISUALIZE_SPROUTS_H
