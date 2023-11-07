//
// Created by david on 11/7/23.
//

#ifndef ARCOS_BLOCK_H
#define ARCOS_BLOCK_H

#include "progargs.hpp"

class block {
public:
    // Constructor
    block(double sx, double sy, double sz, double px, double py, double pz){
      this->sx = sx;
      this->sy = sy;
      this->sz = sz;
      this->px = px;
      this->py = py;
      this->pz = pz;
    }

    // Destructor
    ~block() {
        // Liberación de recursos o acciones de limpieza en el destructor
    }

private:
    /*Tamaño de bloques*/
    double sx;
    double sy;
    double sz;
    double px;
    double py;
    double pz;
};


#endif//ARCOS_BLOCK_H
