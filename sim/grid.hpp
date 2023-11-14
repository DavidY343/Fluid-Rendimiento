//
// Created by david on 11/7/23.
//
#ifndef ARCOS_GRID_HPP
#define ARCOS_GRID_HPP

#include "block.hpp"
#include "progargs.hpp"

#include <iostream>

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
      /* TODO:ver si este comentario se puede eliminar o hay q echar el ojo a algo
       * este bucle no se puede dejar asi, pq luego se recorren todos los
      bloques para ver si las particulas estan bien colocadas
      tb hay q ver como gestionar las particulas con coordenadas fuera de la malla
      */
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

    void simular() {
      // esto lo use para hacer pruebas pero se puede quitar
      //localizar_particulas();

      // 4.3.1
      reposicionar_particulas();

      // 4.3.2
      calcular_aceleraciones();

      // 4.3.3 y 4.3.4 y 4.3.5
      colisiones_particulas();
    }

    //TODO: pendiente eliminar
    //void localizar_particulas();
    void reposicionar_particulas();
    void recolocar_particula(const particula &part);

    void calcular_aceleraciones();
    void inicializar_densidades();
    void calcular_densidades();
    void transformar_densidades();
    void _calcular_aceleraciones();

    void colisiones_particulas();
    void almacenar_resultados(std::ofstream & outputFile);
    void bucle_colisiones(int num_bloque, bool lim_inf, int dimension);

    particula acceder_bloque_part(int b, int p) { return bloques[b].getParticulas()[p]; }

    [[nodiscard]] std::vector<int> obtener_contiguos(int n) const {
      std::vector<int> bloques_contiguos;
      std::vector<int> coordenadas = obtener_coordenadas(n);
      for (int i = -1; i < 2; ++i) {
        for (int j = -1; j < 2; ++j) {
          for (int k = -1; k < 2; ++k) {
            if ((i != 0 || j != 0 || k != 0) &&
                (0 <= coordenadas[0] + i && coordenadas[0] + i < nx && 0 <= coordenadas[1] + j &&
                 coordenadas[1] + j < ny && 0 <= coordenadas[1] + k && coordenadas[1] + k < nz)) {
              bloques_contiguos.push_back(
                  obtener_indice(coordenadas[0] + i, coordenadas[1] + j, coordenadas[1] + k));
            }
          }
        }
      }
      return bloques_contiguos;
    }

    [[nodiscard]] std::vector<int> obtener_coordenadas(int n) const {
      std::vector<int> coordenadas;
      coordenadas.push_back(n / (nz * ny));
      coordenadas.push_back((n % (nz * ny) / nz));
      coordenadas.push_back((n % (nz * ny) % nz));
      return coordenadas;
    }

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
void simulate(std::ifstream const & inputFile, int max_iteraciones, std::ofstream & outputFile);
#endif  // ARCOS_GRID_HPP
