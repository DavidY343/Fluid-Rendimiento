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
  inicializar_densidades();

  // calcular densidades
  calcular_densidades();

  // transformar densidades
  transformar_densidades();

  // calcular aceleraciones
  _calcular_aceleraciones();
}

void grid::inicializar_densidades(){
  for (int indice_bloque = 0; indice_bloque < nx * ny * nz; indice_bloque++) {
    for (auto & particula : bloques[indice_bloque].particulas) {
      particula.inicializar_densidad_aceleracion();
    }
  }
}

void grid::calcular_densidades() {
  for (int indice_b = 0; indice_b < nx * ny * nz; indice_b++) {
    for (unsigned long pi = 0; pi < bloques[indice_b].particulas.size(); pi++) {
      for (unsigned long pj = 0; pj < bloques[indice_b].particulas.size();
           pj++) {  // reducir las iteraciones de este bucle
        if (pi > pj) {
          bloques[indice_b].particulas[pi].interactuar_densidad(bloques[indice_b].particulas[pj], h, true);
        }
      }
    }
  }

  for (int b = 0; b < nx * ny * nz; b++) {
    for (unsigned long pi = 0; pi < bloques[b].particulas.size(); pi++) {
      std::vector<int> const bloques_contiguos = obtener_contiguos(b);
      for (int const bloques_contiguo : bloques_contiguos) {
        for (unsigned long pj = 0; pj < bloques[bloques_contiguo].particulas.size(); pj++) {
          bloques[b].particulas[pi].interactuar_densidad(bloques[bloques_contiguo].particulas[pj], h, false);
        }
      }
    }
  }
}

void grid::transformar_densidades(){
  for (int indice_bloque = 0; indice_bloque < nx * ny * nz; indice_bloque++) {
    for (auto & particula : bloques[indice_bloque].particulas) {
      particula.transformar_densidad(h);
    }
  }
}

void grid::_calcular_aceleraciones(){
  for (int indice_bloque = 0; indice_bloque < nx * ny * nz; indice_bloque++) {
    for (unsigned long pi = 0; pi < bloques[indice_bloque].particulas.size(); pi++) {
      for (unsigned long pj = 0; pj < bloques[indice_bloque].particulas.size();
           pj++) {  // reducir las iteraciones de este bucle
        if (pi > pj) {
          bloques[indice_bloque].particulas[pi].interactuar_aceleracion(bloques[indice_bloque].particulas[pj], h, m,
                                                            true);
        }
      }
    }
  }

  // bloques contiguos
  for (int indice_bloque = 0; indice_bloque < nx * ny * nz; indice_bloque++) {
    for (unsigned long pi = 0; pi < bloques[indice_bloque].particulas.size(); pi++) {
      std::vector<int> const bloques_contiguos = obtener_contiguos(indice_bloque);
      for (unsigned long b2 = 0; b2 < bloques_contiguos.size(); b2++) {
        for (unsigned long pj = 0; pj < bloques[b2].particulas.size(); pj++) {
          bloques[indice_bloque].particulas[pi].interactuar_aceleracion(bloques[b2].particulas[pj], h, m,
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
    }
    if (coordenadas[1] == 0 || coordenadas[1] == getny() - 1) {
      bucle_colisiones(i, coordenadas[1] == 0, 1);
    }
    if (coordenadas[2] == 0 || coordenadas[2] == getnz() - 1) {
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
    cout << "****************************************************\n";
    cout << "iniciando  iteracion " << iteracion << "\n";

    malla.simular();

    cout << "finalizada iteracion " << iteracion << "\n";
    cout << "----------------------------------------------------\n";
  }
}
void escribir_parametros_generales(std::ofstream & outputFile, const std::string& filename) {
  std::ifstream inputFile(filename,
                          std::ios::binary);
  inputFile.seekg(0, std::ios::beg);
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
template <typename T>
void mySwap(T& a, T& b) {
  T temp = a;
  a = b;
  b = temp;
}

template <typename T>
void mySort(std::vector<T>& arr) {
  for (size_t i = 0; i < arr.size() - 1; ++i) {
    for (size_t j = 0; j < arr.size() - 1 - i; ++j) {
        if (arr[j + 1].getid() < arr[j].getid()) {
        mySwap(arr[j], arr[j + 1]);
        }
    }
  }
}
void grid::almacenar_resultados(std::ofstream & outputFile, const std::string& filename) {
  // Escribir los parámetros generales
  escribir_parametros_generales(outputFile, filename);

  std::vector<particula> particulas;
  for (int i = 0; i < getnx() * getny() * getnz(); i++) {
    for (auto & particula : bloques[i].particulas) {
        particulas.push_back(particula);
    }
  }
  mySort(particulas);
  // Escribir los datos de las partículas
  for (const auto & particula : particulas){
    escribir_datos_particulas(outputFile, particula);
  }
}
// No me deja usar mergesort pq me dice que es recursivo xD
/*
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
*/