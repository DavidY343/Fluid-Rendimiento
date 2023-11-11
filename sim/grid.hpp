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
            // cout<<"creando bloque"<<i<<j<<k<<" con indice:"<<obtener_indice(i,j,k)<<"
            // pz="<<constantes::bmin_const[2] + k*sz<<endl;
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
      for (int i = 0; i < 4800; i++) {
        /*
        cout << "creando particula: (" << particulas[i].getpx() << ","
             << particulas[i].getpy() << "," << particulas[i].getpz()
             << ") id=" << particulas[i].getid() << endl;
             */
        // descarta todas las particulas q estan  fuera de la malla (igual luego lo quitamos)
        // if(constantes::bmin_const[0]<=particulas[i].getpx()&&particulas[i].getpx()<constantes::bmax_const[0]&&constantes::bmin_const[1]<=particulas[i].getpy()&&particulas[i].getpy()<constantes::bmax_const[1]&&constantes::bmin_const[2]<=particulas[i].getpz()&&particulas[i].getpz()<constantes::bmax_const[2])
        // {
        recolocar_particula(particulas[i]);
        //}
      }
      cout << "Se inicializo malla con caracteristicas:\n  nx=" << nx << "  ny=" << ny
           << "  nz=" << nz << "\n  sx=" << sx << "  sy=" << sy << "  sz=" << sz << endl;

      int suma = 0;
      for (int b = 0; b < nx * ny * nz; b++) {
        suma += bloques[b].particulas.size();
        // cout<<"bloque:"<<b<<"  con:"<<bloques[b].particulas.size()<<endl;
      }
      cout << "Se encontraron " << suma << " particulas en la malla" << endl;
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

    void localizar_particulas() {
      unsigned long suma = 0;
      for (int b_dat = 0; b_dat < nx * ny * nz; b_dat++) {
        suma += bloques[b_dat].particulas.size();
      }
      using namespace std;
      cout << "Se encontraron " << suma << " particulas en la malla" << endl;
    }

    void reposicionar_particulas() {
      using namespace std;
      cout << "calculando particulas a reposicionar..." << endl;
      std::vector<particula> particulas_a_reposicionar;
      for (int b = 0; b < nx * ny * nz; b++) {
        std::vector<particula> n_part = bloques[b].devolver_particulas();
        if (!n_part.empty()) {
          particulas_a_reposicionar.insert(particulas_a_reposicionar.end(), n_part.begin(),
                                           n_part.end());
        }
      }
      cout << "se reposicionaran: " << particulas_a_reposicionar.size() << " particulas" << endl;
      for (int p = particulas_a_reposicionar.size() - 1; 0 <= p; p--) {
        recolocar_particula(particulas_a_reposicionar[p]);
      }
      cout << "particulas reposicionadas" << endl;
    }

    void recolocar_particula(particula part) {  // esta funcion se puede simplificar
      using namespace std;
      int i = static_cast<int>((part.getpx() - constantes::bmin_const[0]) / sx);
      int j = static_cast<int>((part.getpy() - constantes::bmin_const[1]) / sy);
      int k = static_cast<int>((part.getpz() - constantes::bmin_const[2]) / sz);
      int b = obtener_indice(i, j, k);
      // cout<<i<<","<<j<<","<<k<<","<<b<<"  particula id="<<part.getid()<<endl<<"
      // itam:"<<bloques[b].particulas.size();
      bloques[b].anhadir_particulas(part);
      // cout<<"ntam:"<<bloques[b].particulas.size()<<endl;
      // cout<<"particula: ("<<part.getpx()<<","<<part.getpy()<<","<<part.getpz()<<")\n";
      // cout<<" en bloque:"<<i<<","<<j<<","<<k<<" "<<b<<endl;
      // cout<<bloques[b].getpz()<<endl;
    }

    void calcular_aceleraciones() {
      /* seguimiento de una particula, para hacer pruebas
      for(int b=0; b<nx*ny*nz;b++) {
        for(unsigned long pi=0; pi<bloques[b].particulas.size();pi++){
          if(bloques[b].particulas[pi].getid() == 3787) { //3787 tb sirve para las pruebas
          bloques[b].particulas[pi].imprimir_datos();
          }
        }
      }
       */

      // inicializar
      for (int b = 0; b < nx * ny * nz; b++) {
        for (unsigned long p = 0; p < bloques[b].particulas.size(); p++) {
          bloques[b].particulas[p].inicializar_densidad_aceleracion();
        }
      }

      // calcular densidades, iteraciones con particulas del mismo bloque
      for (int b = 0; b < nx * ny * nz; b++) {
        for (unsigned long pi = 0; pi < bloques[b].particulas.size(); pi++) {
          for (unsigned long pj = 0; pj < bloques[b].particulas.size();
               pj++) {  // reducir las iteraciones de este bucle
            if (pi > pj) {
              bloques[b].particulas[pi].interactuar_densidad(bloques[b].particulas[pj], h, true);
            }
          }
        }
      }

      // bloques contiguos
      for (int b = 0; b < nx * ny * nz; b++) {
        for (unsigned long pi = 0; pi < bloques[b].particulas.size(); pi++) {
          std::vector<int> bloques_contiguos = obtener_contiguos(b);
          for (unsigned long b2 = 0; b2 < bloques_contiguos.size(); b2++) {
            for (unsigned long pj = 0; pj < bloques[b2].particulas.size(); pj++) {
              bloques[b].particulas[pi].interactuar_densidad(bloques[b2].particulas[pj], h, false);
            }
          }
        }
      }

      // seguimiento de una particula, para hacer pruebas
      /*
      for(int b=0; b<1;b++) {
        for(unsigned long pi=0; pi<bloques[b].particulas.size();pi++){
          //if(bloques[b].particulas[pi].getid() == 0) { //3787 tb sirve para las pruebas
            bloques[b].particulas[pi].imprimir_datos();
          //}
        }
      }
       */

      // transformar densidades
      for (int b = 0; b < nx * ny * nz; b++) {
        for (unsigned long p = 0; p < bloques[b].particulas.size(); p++) {
          bloques[b].particulas[p].transformar_densidad(h);
        }
      }

      // calcular aceleraciones, iteraciones con particulas del mismo bloque
      for (int b = 0; b < nx * ny * nz; b++) {
        for (unsigned long pi = 0; pi < bloques[b].particulas.size(); pi++) {
          for (unsigned long pj = 0; pj < bloques[b].particulas.size();
               pj++) {  // reducir las iteraciones de este bucle
            if (pi > pj) {
              bloques[b].particulas[pi].interactuar_aceleracion(bloques[b].particulas[pj], h, m,
                                                                true);
            }
          }
        }
      }

      // bloques contiguos
      for (int b = 0; b < nx * ny * nz; b++) {
        for (unsigned long pi = 0; pi < bloques[b].particulas.size(); pi++) {
          std::vector<int> bloques_contiguos = obtener_contiguos(b);
          for (unsigned long b2 = 0; b2 < bloques_contiguos.size(); b2++) {
            for (unsigned long pj = 0; pj < bloques[b2].particulas.size(); pj++) {
              bloques[b].particulas[pi].interactuar_aceleracion(bloques[b2].particulas[pj], h, m,
                                                                false);
            }
          }
        }
      }
    }

    void colisiones_particulas();

    void bucle_colisiones(int num_bloque, bool lim_inf, int dimension);

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

void init_simulate(int max_iteraciones, grid & malla);
grid init_params(std::ifstream const & inputFile);
#endif  // ARCOS_GRID_HPP
