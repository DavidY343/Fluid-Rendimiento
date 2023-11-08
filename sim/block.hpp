//
// Created by david on 11/7/23.
//

#ifndef ARCOS_BLOCK_H
#define ARCOS_BLOCK_H

#include <vector>
#include "progargs.hpp"
#include <iostream>

class block {
public:

    // Constructor
    block(double px, double py, double pz, double sx, double sy, double sz){
      this->sx = sx;
      this->sy = sy;
      this->sz = sz;
      this->px = px;
      this->py = py;
      this->pz = pz;
    }
    std::vector<particula> particulas;

    // Destructor
    ~block() {
        // Liberación de recursos o acciones de limpieza en el destructor
    }

    std::vector<particula> devolver_particulas(){
      std::vector<particula> devolver;
      for(int p = particulas.size()-1; 0<=p; p--){
        if (!p_bloque(particulas[p])){

             devolver.push_back(particulas[p]);
             //particulas.erase(particulas.begin() + p);
             this->particulas = eliminar(particulas, p);

        }
      }
      return devolver;
    }

    //pruebas
    void imprimir_vector(std::vector<particula> v){
      using namespace std;
      for(unsigned long i = 0; i<v.size(); i++) {
        cout << "v["<<i<<"]: "<<v[i].getid()<<endl;
      }
    }

    //espero poder elinar esta funcion
    std::vector<particula> eliminar(std::vector<particula> v, int e){
      std::vector<particula> eliminado;
      for(unsigned long i = 0; i<v.size(); i++) {
        if(e != static_cast<int>(i))
          eliminado.push_back(v[i]);
      }
      return eliminado;
    }

    void anhadir_particulas(particula part){
        particulas.push_back(part);
    }

    [[nodiscard]] bool p_bloque(particula p) const {
      return px<=p.getpx()&&p.getpx()<px+sx&&py<=p.getpy()&&p.getpy()<py+sy&&pz<=p.getpz()&&p.getpz()<pz+sz;
    }

    [[nodiscard]] double getpz() const{
      return pz;
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
