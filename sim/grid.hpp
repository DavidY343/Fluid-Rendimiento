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

      using namespace std;

      sx = (constantes::bmax_const[0] - constantes::bmin_const[0]) / nx;
      sy = (constantes::bmax_const[1] - constantes::bmin_const[1]) / ny;
      sz = (constantes::bmax_const[2] - constantes::bmin_const[2]) / nz;
      for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
          for (int k = 0; k < nz; k++) {
            //cout<<"creando bloque"<<i<<j<<k<<" con indice:"<<obtener_indice(i,j,k)<<" pz="<<constantes::bmin_const[2] + k*sz<<endl;
            bloques.emplace_back(constantes::bmin_const[0] + i*sx,constantes::bmin_const[1] + j*sy, constantes::bmin_const[2] + k*sz, sx, sy, sz);
            if(obtener_indice(i,j,k) >= 23){
            }
          }
        }
      }

      /* este bucle no se puede dejar asi, pq luego se recorren todos los
      bloques para ver si las particulas estan bien colocadas
      tb hay q ver como gestionar las particulas con coordenadas fuera de la malla */
      for(int i = 0; i<4800;i++){
        //cout<<"creando particula: ("<<particulas[i].getpx()<<","<<particulas[i].getpy()<<","<<particulas[i].getpz()<<")"<<endl;
        // descarta todas las particulas q estan  fuera de la malla (igual luego lo quitamos)
        if(constantes::bmin_const[0]<=particulas[i].getpx()&&particulas[i].getpx()<constantes::bmax_const[0]&&constantes::bmin_const[1]<=particulas[i].getpy()&&particulas[i].getpy()<constantes::bmax_const[1]&&constantes::bmin_const[2]<=particulas[i].getpz()&&particulas[i].getpz()<constantes::bmax_const[2]) {
          recolocar_particula(particulas[i]);
        }
      }
      cout<<"Se inicializo malla con caracteristicas:\n  nx="<<nx<<"  ny="<<ny<<"  nz="<<nz<<"\n  sx="<<sx<<"  sy="<<sy<<"  sz="<<sz<<endl;
    }

    // Destructor
    ~grid() {
      // Liberación de recursos o acciones de limpieza en el destructor
    }

    void simular(){
      //4.3.1
      reposicionar_particulas();

      //4.3.2
      calcular_aceleraciones();

    }

    void reposicionar_particulas(){
      using namespace std;
      cout<<"calculando particulas a reposicionar"<<endl;
      std::vector<particula> particulas_a_reposicionar;
      for(int b=0; b<nx*ny*nz;b++) {
        std::vector<particula> n_part = bloques[b].devolver_particulas();
        if(!n_part.empty()){
          particulas_a_reposicionar.insert(particulas_a_reposicionar.end(), n_part.begin(), n_part.end());
        }
      }
      cout<<"se reposicionaran: "<<particulas_a_reposicionar.size()<<" particulas"<<endl;
      for(int p = particulas_a_reposicionar.size()-1; 0<=p; p--){
        recolocar_particula(particulas_a_reposicionar[p]);
      }
      cout<<"particulas reposicionadas"<<endl;
    }

    void recolocar_particula(particula part){ //esta funcion se puede simplificar
      using namespace std;
      int i = static_cast<int>((part.getpx()-constantes::bmin_const[0])/sx);
      int j = static_cast<int>((part.getpy()-constantes::bmin_const[1])/sy);
      int k = static_cast<int>((part.getpz()-constantes::bmin_const[2])/sz);
      int b = obtener_indice(i, j, k);
      bloques[b].anhadir_particulas(part);
      //cout<<"particula: ("<<part.getpx()<<","<<part.getpy()<<","<<part.getpz()<<")";
      //cout<<" en bloque:"<<i<<","<<j<<","<<k<<" "<<b<<endl;
      //cout<<bloques[b].getpz()<<endl;
    }

    void calcular_aceleraciones(){
      //inicializar
      for(int b=0; b<nx*ny*nz;b++) {
        for(unsigned long p=0; p<bloques[b].particulas.size();p++){
          bloques[b].particulas[p].inicializar_densidad_aceleracion();
        }
      }

      //calcular densidades, iteraciones con particulas del mismo bloque
      for(int b=0; b<nx*ny*nz;b++) {
        for(unsigned long pi=0; pi<bloques[b].particulas.size();pi++){
          for(unsigned long pj=0; pj<bloques[b].particulas.size();pj++) //reducir las iteraciones de este bucle
            if(pi>pj) {
              bloques[b].particulas[pi].interactuar_densidad(bloques[b].particulas[pj], h); //esto no aplica a pj la interaccion
            }
        }
      }
      //bloques contiguos


      //transformar densidades
      for(int b=0; b<nx*ny*nz;b++) {
        for(unsigned long p=0; p<bloques[b].particulas.size();p++){
          bloques[b].particulas[p].transformar_densidad(h);
        }
      }

      //calcular aceleraciones, iteraciones con particulas del mismo bloque
      for(int b=0; b<nx*ny*nz;b++) {
        for(unsigned long pi=0; pi<bloques[b].particulas.size();pi++){
          for(unsigned long pj=0; pj<bloques[b].particulas.size();pj++) //reducir las iteraciones de este bucle
            if(pi>pj) {
              bloques[b].particulas[pi].interactuar_aceleracion(bloques[b].particulas[pj], h, m); //esto no aplica a pj la interaccion
            }
        }
      }
    }

    [[nodiscard]] int obtener_indice(int i, int j, int k) const{
      return nz*ny*i + nz*j + k;
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
    double sx;
    double sy;
    double sz;
    std::vector<block> bloques;
};


#endif//ARCOS_GRID_H
