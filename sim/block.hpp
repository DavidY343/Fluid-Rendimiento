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
    block(std::vector<double> const &position, std::vector<double> const &size) : sx(size[0]), sy(size[1]),sz(size[2]), px(position[0]), py(position[1]),pz(position[2]){
    }
    // Constructor de copia
    block(const block& other) = default;

    // Operador de asignación de copia
    block& operator=(const block& other) {
      if (this != &other) { // Check for self-assignment
        sx = other.sx;
        sy = other.sy;
        sz = other.sz;
        px = other.px;
        py = other.py;
        pz = other.pz;
        particulas = other.particulas;
      }
      return *this;
    }

    // Constructor de movimiento
    block(block && other) = default;

    // Destructor
    ~block() = default;

    // Operador de asignación de movimiento
    block & operator=(block && other) noexcept {
      if (this != &other) { }
      return *this;
    }
    std::vector<particula> devolver_particulas(){
      std::vector<particula> devolver;
      for(int p_dat = particulas.size()-1; 0<=p_dat; p_dat--){
        if (!p_bloque(particulas[p_dat])){

          devolver.push_back(particulas[p_dat]);
          //particulas.erase(particulas.begin() + p);
          this->particulas = eliminar(particulas, p_dat);

        }
      }
      return devolver;
    }


    //espero poder elinar esta funcion
    static std::vector<particula> eliminar(std::vector<particula> v, int e){
      std::vector<particula> eliminado;
      for(unsigned long i = 0; i<v.size(); i++) {
        if(e != static_cast<int>(i)){
          eliminado.push_back(v[i]);}
      }
      return eliminado;
    }
    void anhadir_particulas(const particula &part){
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

public:
  [[nodiscard]] double getSx() const { return sx; }

  [[nodiscard]] double getSy() const { return sy; }

  [[nodiscard]] double getSz() const { return sz; }

  [[nodiscard]] double getPx() const { return px; }

  [[nodiscard]] double getPy() const { return py; }

  [[nodiscard]] double getPz() const { return pz; }

  [[nodiscard]] std::vector<particula> & getParticulas() { return particulas; }

  void setSx(double new_sx) {sx = new_sx; }

  void setSy(double new_sy) {sy = new_sy; }

  void setSz(double new_sz) {sz = new_sz; }

  void setPx(double new_px) {px = new_px; }

  void setPy(double new_py) {py = new_py; }

  void setPz(double new_pz) {pz = new_pz; }

  void setParticulas(std::vector<particula> & new_particulas) {particulas = new_particulas; }

private:
  double sy;
  double sz;
  double px;
  double py;
  double pz;
  std::vector<particula> particulas;
};


#endif//ARCOS_BLOCK_H
