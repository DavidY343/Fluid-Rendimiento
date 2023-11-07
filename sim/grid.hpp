//
// Created by david on 11/7/23.
//

#ifndef ARCOS_GRID_HPP
#define ARCOS_GRID_HPP


class grid {
  public:
    // Constructor
    grid(double ppm, int np, int nx, int ny, int nz) : ppm(ppm), np(np), nx(nx), ny(ny), nz(nz) { }

    // Destructor
    ~grid() {
      // Liberación de recursos o acciones de limpieza en el destructor
    }

    /*Getters*/
    [[nodiscard]] int getnx() const { return nx; }

    [[nodiscard]] int getny() const { return ny; }

    [[nodiscard]] int getnz() const { return nz; }

  private:
    double ppm;
    int np;
    /*Número de bloques*/
    int nx;
    int ny;
    int nz;
};


#endif//ARCOS_GRID_H
