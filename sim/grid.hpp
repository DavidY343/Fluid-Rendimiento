//
// Created by david on 11/7/23.
//
#ifndef ARCOS_GRID_HPP
#define ARCOS_GRID_HPP

#include <iostream>
#include "progargs.hpp"
#include "block.hpp"

class grid {
  public:
    // Constructor
    grid(std::vector<double> cabeceras, std::vector<particula> particulas) {
      m = cabeceras[0];
      h = cabeceras[1];
      nx = static_cast<int>((constantes::bmax_const[0] - constantes::bmin_const[0]) / h);
      ny = static_cast<int>((constantes::bmax_const[1] - constantes::bmin_const[1]) / h);
      nz = static_cast<int>((constantes::bmax_const[2] - constantes::bmin_const[2]) / h);

      double sx = (constantes::bmax_const[0] - constantes::bmin_const[0]) / nx;
      double sy = (constantes::bmax_const[1] - constantes::bmin_const[1]) / ny;
      double sz = (constantes::bmax_const[2] - constantes::bmin_const[2]) / nz;
      for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
          for (int k = 0; k < nz; k++) {
            bloques.emplace_back(nx, ny, nz, i * sx, j * sy, k * sz);
          }
        }
      }
      using namespace std;
      cout<<particulas[0].getpx();
    }

    // Destructor
    ~grid() {
      // Liberación de recursos o acciones de limpieza en el destructor
    }

    /*Getters*/
    [[nodiscard]] int getnx() const { return nx; }

    [[nodiscard]] int getny() const { return ny; }

    [[nodiscard]] int getnz() const { return nz; }

  private:
    double m;
    double h;
    /*Número de bloques*/
    int nx;
    int ny;
    int nz;
    std::vector<block> bloques;
};


#endif//ARCOS_GRID_H
