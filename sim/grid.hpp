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
    grid(std::vector<double> cabeceras, std::vector<particula> particulas) {
      m  = cabeceras[0];
      h  = cabeceras[1];
      nx = static_cast<int>((constantes::bmax_const[0] - constantes::bmin_const[0]) / h);
      ny = static_cast<int>((constantes::bmax_const[1] - constantes::bmin_const[1]) / h);
      nz = static_cast<int>((constantes::bmax_const[2] - constantes::bmin_const[2]) / h);

      using namespace std;

      sx = (constantes::bmax_const[0] - constantes::bmin_const[0]) / nx;
      sy = (constantes::bmax_const[1] - constantes::bmin_const[1]) / ny;
      sz = (constantes::bmax_const[2] - constantes::bmin_const[2]) / nz;
      for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
          for (int k = 0; k < nz; k++) {
            bloques.emplace_back(constantes::bmin_const[0] + i * sx,
                                 constantes::bmin_const[1] + j * sy,
                                 constantes::bmin_const[2] + k * sz, sx, sy, sz);
          }
        }
      }

      /* este bucle no se puede dejar asi, pq luego se recorren todos los
      bloques para ver si las particulas estan bien colocadas
      tb hay q ver como gestionar las particulas con coordenadas fuera de la malla
    */
      for (const auto & particula : particulas) {
        // descarta todas las particulas q estan  fuera de la malla (igual luego lo quitamos)
        // if(constantes::bmin_const[0]<=particulas[i].getpx()&&particulas[i].getpx()<constantes::bmax_const[0]&&constantes::bmin_const[1]<=particulas[i].getpy()&&particulas[i].getpy()<constantes::bmax_const[1]&&constantes::bmin_const[2]<=particulas[i].getpz()&&particulas[i].getpz()<constantes::bmax_const[2])
        // {
        recolocar_particula(particula);
        //}
      }
    }

    // Destructor
    ~grid() {
      // Liberación de recursos o acciones de limpieza en el destructor
    }

    void simular() {
      // esto lo use para hacer pruebas pero se puede quitar
      localizar_particulas();

      // 4.3.1
      reposicionar_particulas();

      // 4.3.2
      calcular_aceleraciones();

      // 4.3.3 y 4.3.4 y 4.3.5
      colisiones_particulas();
    }

    void localizar_particulas();

    void reposicionar_particulas();

    void recolocar_particula(particula part);

    void calcular_aceleraciones();

    void colisiones_particulas();
    void almacenar_resultados(std::ofstream & outputFile, std::ifstream const & inputFile);
    void bucle_colisiones(int num_bloque, bool lim_inf, int dimension);

    int acceder_bloque_part(int b, int p){
      return bloques[b].particulas[p].getid();
    }

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
    convertirDatos(particula const & particula) ;
        void init_simulate(int max_iteraciones, grid & malla);
grid init_params(std::ifstream const & inputFile);
void bubbleSort(std::vector<particula> &particles);
#endif  // ARCOS_GRID_HPP
