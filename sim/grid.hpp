//
// Created by david on 11/7/23.
//
#ifndef ARCOS_GRID_HPP
#define ARCOS_GRID_HPP

#include "block.hpp"
#include "progargs.hpp"

#include <iostream>
#include <sstream>

class grid {
  public:
    // Constructor
    grid(std::vector<double> cabeceras, std::vector<particula> const & particulas)
      : m(cabeceras[0]), h(cabeceras[1]),
        nx(static_cast<int>((constantes::bmax_const[0] - constantes::bmin_const[0]) / h)),
        ny(static_cast<int>((constantes::bmax_const[1] - constantes::bmin_const[1]) / h)),
        nz(static_cast<int>((constantes::bmax_const[2] - constantes::bmin_const[2]) / h)),
        sx((constantes::bmax_const[0] - constantes::bmin_const[0]) / nx),
        sy((constantes::bmax_const[1] - constantes::bmin_const[1]) / ny),
        sz((constantes::bmax_const[2] - constantes::bmin_const[2]) / nz) {
      using namespace std;
      for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
          for (int k = 0; k < nz; k++) {
            std::vector<double> const position = {constantes::bmin_const[0] + i * sx,
                                                  constantes::bmin_const[1] + j * sy,
                                                  constantes::bmin_const[2] + k * sz};
            std::vector<double> const size = {sx, sy, sz};
            bloques.emplace_back(position,size);
          }
        }
      }
      for (auto const & particula : particulas) {
        recolocar_particula(particula);
      }
    }

    grid(grid const & other) = default;

    // Copy assignment operator
    grid & operator=(grid const & other) {
      if (this != &other) {
        m       = other.m;
        h       = other.h;
        nx      = other.nx;
        ny      = other.ny;
        nz      = other.nz;
        sx      = other.sx;
        sy      = other.sy;
        sz      = other.sz;
        bloques = other.bloques;
      }
      return *this;
    }

    // Move constructor
    grid(grid && other) noexcept
      : m(other.m), h(other.h), nx(other.nx), ny(other.ny), nz(other.nz), sx(other.sx),
        sy(other.sy), sz(other.sz), bloques(std::move(other.bloques)) { }

    // Move assignment operator
    grid & operator=(grid && other) noexcept {
      if (this != &other) {
        m       = other.m;
        h       = other.h;
        nx      = other.nx;
        ny      = other.ny;
        nz      = other.nz;
        sx      = other.sx;
        sy      = other.sy;
        sz      = other.sz;
        bloques = std::move(other.bloques);
      }
      return *this;
    }

    // Destructor
    ~grid() = default;

    void simular();

    void reposicionar_particulas();
    void recolocar_particula(const particula &part);

    void inicializar_densidades();
    void calcular_densidades();
    void transformar_densidades();
    void transferir_aceleraciones();

    void colisiones_particulas();
    void movimiento_particulas();
    void bucle_limites(int num_bloque, bool lim_inf, int dimension);
    void almacenar_resultados(std::ofstream & outputFile, double ppm, std::vector<particula> const & part);
    void bucle_colisiones(int num_bloque, bool lim_inf, int dimension);

    particula acceder_bloque_part(int b, int p) { return bloques[b].getParticulas()[p]; }

    [[nodiscard]] std::vector<int> obtener_contiguos(int n) const;
    [[nodiscard]] std::vector<int> obtener_coordenadas(int n) const;
    [[nodiscard]] int obtener_indice(int i, int j, int k) const { return nz * ny * i + nz * j + k; }

    /*Getters*/
    [[nodiscard]] int getnx() const { return nx; }

    [[nodiscard]] int getny() const { return ny; }

    [[nodiscard]] int getnz() const { return nz; }

    [[nodiscard]] double getsx() const { return sx; }

    [[nodiscard]] double getsy() const { return sy; }

    [[nodiscard]] double getsz() const { return sz; }

  private:
    double m;
    double h;
    /*Número de bloques*/
    int nx;
    int ny;
    int nz;
    double sx;
    double sy;
    double sz;
    std::vector<block> bloques;
};

void escribir_datos_particulas(std::ofstream & outputFile, particula const & particula);
std::tuple<float, float, float, float, float, float, float, float, float>
    convertirDatos(particula const & particula);
void init_simulation(std::ifstream const & inputFile, int max_iteraciones, std::ofstream & outputFile);
#endif  // ARCOS_GRID_HPP
