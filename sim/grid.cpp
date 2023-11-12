//
// Created by david on 11/7/23.
//
//
// Created by david on 10/3/23.
//

#include "grid.hpp"
#include "progargs.hpp"
#include "progargs.cpp"
#include <iostream>
#include <tuple>

void grid::localizar_particulas(){
  unsigned long suma = 0;
  for (int b_dat = 0; b_dat < nx * ny * nz; b_dat++) {
    suma += bloques[b_dat].particulas.size();
  }
  using namespace std;
  cout << "Se encontraron " << suma << " particulas en la malla\n";
}

void grid::reposicionar_particulas() {
  using namespace std;
  cout << "calculando particulas a reposicionar...\n";
  std::vector<particula> particulas_a_reposicionar;
  for (int indice_bloque = 0; indice_bloque < nx * ny * nz; indice_bloque++) {
    std::vector<particula> n_part = bloques[indice_bloque].devolver_particulas();
    if (!n_part.empty()) {
      particulas_a_reposicionar.insert(particulas_a_reposicionar.end(), n_part.begin(),
                                       n_part.end());
    }
  }
  cout << "se reposicionaran: " << particulas_a_reposicionar.size() << " particulas\n";
  for (int par = static_cast<int>(particulas_a_reposicionar.size()) - 1; 0 <= par; par--) {
    recolocar_particula(particulas_a_reposicionar[par]);
  }
  cout << "particulas reposicionadas\n";
}

void grid::recolocar_particula(particula part) {  // esta funcion se puede simplificar
  using namespace std;
  int indice_i = static_cast<int>((part.getpx() - constantes::bmin_const[0]) / sx);
  if(indice_i<0){
    indice_i = 0;
  } else if (indice_i >= nx){
    indice_i = nx-1;
  }
  int indice_j = static_cast<int>((part.getpy() - constantes::bmin_const[1]) / sy);
  if(indice_j<0){
    indice_j = 0;
  } else if (indice_j >= ny){
    indice_j = ny-1;
  }
  int indice_k = static_cast<int>((part.getpz() - constantes::bmin_const[2]) / sz);
  if(indice_k<0){
    indice_k = 0;
  } else if (indice_k >= nz){
    indice_k = nz-1;
  }
  int const indice_bloque = obtener_indice(indice_i, indice_j, indice_k);
  //cout<<indice_i<<","<<indice_j<<","<<indice_k<<","<<indice_bloque<<"  particula id="<<part.getid()<<endl;
  bloques[indice_bloque].anhadir_particulas(part);
}

void grid::calcular_aceleraciones() {
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

void grid::colisiones_particulas() {
  std::vector<int> coordenadas;
  for (int i = 0; i < getnx() * getny() * getnz(); i++) {
    coordenadas = obtener_coordenadas(i);
    if (coordenadas[0] == 0 || coordenadas[0] == getnx() - 1) {
      bucle_colisiones(i, coordenadas[0] == 0, 0);
    } else if (coordenadas[1] == 0 || coordenadas[1] == getny() - 1) {
      bucle_colisiones(i, coordenadas[1] == 0, 1);
    } else if (coordenadas[2] == 0 || coordenadas[2] == getnz() - 1) {
      bucle_colisiones(i, coordenadas[1] == 0, 2);
    }
  }
}

void grid::bucle_colisiones(int num_bloque, bool lim_inf, int dimension) {
  for (auto & particula : bloques[num_bloque].particulas) {
      if (dimension == 0) {
        particula.colisionLimiteEjeX(lim_inf);
        // 4.3.4
        particula.actualizarMovimiento();
        // 4.3.5
        particula.limiteRecintox(lim_inf);
      }else if (dimension == 1) {
        particula.colisionLimiteEjeY(lim_inf);
        // 4.3.4
        particula.actualizarMovimiento();
        // 4.3.5
        particula.limiteRecintoy(lim_inf);
      }else if (dimension == 2) {
        particula.colisionLimiteEjeZ(lim_inf);
        // 4.3.4
        particula.actualizarMovimiento();
        // 4.3.5
        particula.limiteRecintoz(lim_inf);
      }
  }
}

grid init_params(std::ifstream const & inputFile) {
  auto ppm = static_cast<double>(
      read_binary_value<float>((std::istream &) inputFile));
  double const m_dat = constantes::p_const / std::pow(ppm, 3.0);
  double const h_dat = constantes::r_const / ppm;
  std::vector<double> vtor                = {m_dat, h_dat};
  std::vector<particula> const particulas = crear_particulas(inputFile);
  using namespace std;
  grid const malla(vtor, particulas);
  cout << "Number of particles: " << particulas.size() << "\n";
  cout << "Particles per meter: " << ppm << "\n";
  cout << "Smoothing length: " << vtor[1] << "\n";
  cout << "Particle mass: " << vtor[0] << "\n";
  return malla;
}

void init_simulate(int const max_iteraciones, grid & malla) {
  using namespace std;
  cout << "Grid size: " << malla.getnx() << " x " << malla.getny()  << " x " << malla.getnz() << "\n";
  cout << "Number of blocks: " << malla.getnx() * malla.getny() * malla.getnz() << "\n";
  cout << "Block size: " << malla.getsx() << " x " << malla.getsy()  << " x " << malla.getsz() << "\n";
  for (int iteracion = 1; iteracion <= max_iteraciones; iteracion++) {
    cout << "****************************************************" << endl;
    cout << "iniciando  iteracion " << iteracion << endl;

    malla.simular();

    cout << "finalizada iteracion " << iteracion << endl;
    cout << "----------------------------------------------------" << endl;
  }
}
void escribir_parametros_generales(std::ofstream & outputFile, std::ifstream const & inputFile) {
  auto ppm = (read_binary_value<float>((std::istream &) inputFile));
  auto n_particulas_int = read_binary_value<int>((std::istream &) inputFile);
  outputFile.write(as_buffer(ppm), sizeof(ppm));
  outputFile.write(as_buffer(n_particulas_int), sizeof(n_particulas_int));
}

// Función para convertir los datos de las partículas de doble a simple precisión
std::tuple<float, float, float, float, float, float, float, float, float>
    convertirDatos(particula const & particula) {
  auto px_dat = static_cast<float>(particula.getpx());
  auto py_dat = static_cast<float>(particula.getpy());
  auto pz_dat = static_cast<float>(particula.getpz());
  auto hvx    = static_cast<float>(particula.gethvx());
  auto hvy    = static_cast<float>(particula.gethvy());
  auto hvz    = static_cast<float>(particula.gethvz());
  auto vx_dat = static_cast<float>(particula.getvx());
  auto vy_dat = static_cast<float>(particula.getvy());
  auto vz_dat = static_cast<float>(particula.getvz());

  return std::make_tuple(px_dat, py_dat, pz_dat, hvx, hvy, hvz, vx_dat, vy_dat, vz_dat);
}

void escribir_datos_particulas(std::ofstream & outputFile, particula const & particula) {
  auto [px_dat, py_dat, pz_dat, hvx, hvy, hvz, vx_dat, vy_dat, vz_dat] = convertirDatos(particula);

  outputFile.write(as_buffer(px_dat), sizeof(px_dat));
  outputFile.write(as_buffer(py_dat), sizeof(py_dat));
  outputFile.write(as_buffer(pz_dat), sizeof(pz_dat));
  outputFile.write(as_buffer(hvx), sizeof(hvx));
  outputFile.write(as_buffer(hvy), sizeof(hvy));
  outputFile.write(as_buffer(hvz), sizeof(hvz));
  outputFile.write(as_buffer(vx_dat), sizeof(vx_dat));
  outputFile.write(as_buffer(vy_dat), sizeof(vy_dat));
  outputFile.write(as_buffer(vz_dat), sizeof(vz_dat));
}

/*
void grid::almacenar_resultados(std::ofstream & outputFile, std::ifstream const & inputFile) {
  // Escribir los parámetros generales
  escribir_parametros_generales(outputFile, inputFile);


  for (int i = 0; i < getnx() * getny() * getnz(); i++) {
    for (auto & particula : bloques[i].particulas) {
        escribir_datos_particulas(outputFile, particula);
    }
  }
}*/

void merge(std::vector<particula>& arr, size_t l, size_t m, size_t r) {
  size_t n1 = m - l + 1;
  size_t n2 = r - m;

  // Crear vectores temporales
  std::vector<particula> L(arr.begin() + l, arr.begin() + l + n1);
  std::vector<particula> R(arr.begin() + m + 1, arr.begin() + m + 1 + n2);

  // Índices iniciales de los subvectores
  size_t i = 0, j = 0, k = l;

  // Combinar los subvectores de nuevo en arr
  while (i < n1 && j < n2) {
    if (L[i].getid() <= R[j].getid()) {
        arr[k++] = L[i++];
    } else {
        arr[k++] = R[j++];
    }
  }

  // Copiar los elementos restantes de L (si los hay)
  while (i < n1) {
    arr[k++] = L[i++];
  }

  // Copiar los elementos restantes de R (si los hay)
  while (j < n2) {
    arr[k++] = R[j++];
  }
}

// Función principal de merge sort
void mergeSort(std::vector<particula>& arr, size_t l, size_t r) {
  if (l < r) {
    // Encuentra el punto medio
    size_t m = l + (r - l) / 2;

    // Ordena la primera y la segunda mitad
    mergeSort(arr, l, m);
    mergeSort(arr, m + 1, r);

    // Combina las mitades ordenadas
    merge(arr, l, m, r);
  }
}
bool compareByParticleId(const particula &a, const particula &b) {
  return a.getid() < b.getid();
}

void grid::almacenar_resultados(std::ofstream & outputFile, std::ifstream const & inputFile) {
  // Escribir los parámetros generales
  escribir_parametros_generales(outputFile, inputFile);

  // Escribir los datos de las partículas

  for (int i = 0; i < getnx() * getny() * getnz(); i++) {
    for (auto & particula : bloques[i].particulas) {
        escribir_datos_particulas(outputFile, particula);
    }
  }
}/*
void grid::almacenar_resultados(std::ofstream & outputFile, std::ifstream const & inputFile) {
  // Escribir los parámetros generales
  escribir_parametros_generales(outputFile, inputFile);

  // Escribir los datos de las partículas
  for (int i = 0; i < getnx() * getny() * getnz(); i++) {
    // Ordenar las partículas en el bloque actual
    size_t n = bloques[i].particulas.size();
    for (size_t j = 0; j < n - 1; ++j) {
        for (size_t k = 0; k < n - j - 1; ++k) {
          if (bloques[i].particulas[k].getid() > bloques[i].particulas[k + 1].getid()) {
            // Swap only if indices are within bounds
            if (k + 1 < n) { std::swap(bloques[i].particulas[k], bloques[i].particulas[k + 1]); }
            }
        }
    }
    for (auto & particula : bloques[i].particulas) {
        escribir_datos_particulas(outputFile, particula);
    }
  }
}
*/