//
// Created by david on 11/7/23.
//
//
// Created by david on 10/3/23.
//

#include "grid.hpp"
/*
// Función para colisiones con límites en el eje x
void interaccionesLimitesEjeX(Particula& particula, int cx, double xmin, double xmax) {
  if (cx == 0 || cx == nx - 1) {
    double dx = (cx == 0) ? particula.px - xmin : xmax - particula.px;
    if (dx < 0) {
      particula.px = (cx == 0) ? xmin - dx : xmax + dx;
      particula.vx = -particula.vx;
      particula.hvx = -particula.hvx;
    }
  }
}

// Función para colisiones con límites en el eje y
void interaccionesLimitesEjeY(Particula& particula, int cy, double ymin, double ymax) {
  if (cy == 0 || cy == ny - 1) {
    double dy = (cy == 0) ? particula.py - ymin : ymax - particula.py;
    if (dy < 0) {
      particula.py = (cy == 0) ? ymin - dy : ymax + dy;
      particula.vy = -particula.vy;
      particula.hvy = -particula.hvy;
    }
  }
}

// Función para colisiones con límites en el eje z
void interaccionesLimitesEjeZ(Particula& particula, int cz, double zmin, double zmax) {
  if (cz == 0 || cz == nz - 1) {
    double dz = (cz == 0) ? particula.pz - zmin : zmax - particula.pz;
    if (dz < 0) {
      particula.pz = (cz == 0) ? zmin - dz : zmax + dz;
      particula.vz = -particula.vz;
      particula.hvz = -particula.hvz;
    }
  }
}

// Función general que llama a las tres funciones anteriores
void interaccionesLimitesRecinto(Particula& particula, int cx, int cy, int cz, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
  interaccionesLimitesEjeX(particula, cx, xmin, xmax);
  interaccionesLimitesEjeY(particula, cy, ymin, ymax);
  interaccionesLimitesEjeZ(particula, cz, zmin, zmax);
}

int calculoCx(Particula& particula, int nx, double sx) {
  int cx = static_cast<int>(particula.px / sx);
  if (cx < 0) {
    cx = 0;
  }else if (cx > nx - 1) {
    cx = nx - 1;
  }
  return cx;
}

int calculoCy(Particula& particula, int ny, double sy) {
  int cy = static_cast<int>(particula.py / sy);
  if (cy < 0) {
    cy = 0;
  }else if (cy > ny - 1) {
    cy = ny - 1;
  }
  return cy;
}

int calculoCy(Particula& particula, int nz, double sz) {
  int cz = static_cast<int>(particula.pz / sz);
  if (cz < 0) {
    cz = 0;
  }else if (cz > nz - 1) {
    cz = nz - 1;
  }
  return cz;
}

double colisionLimiteEjeX(Particula& particula, double xmin, int cx, double xmax, double tiempo, double dp_var, int nx, double ax, double ps, double dv) {
  double nuevaX = particula.px + particula.hvx * tiempo;
  double difLimX = 0;
  double const min_value = 0.0000000001;
  double nuevaAx = ax;
    if (cx == 0){
      difLimX = dp_var - (particula.px - xmin);
      if (difLimX > min_value){
        nuevaAx = ax + (ps * difLimX - dv * particula.vx);
      }
    }else if(cx == nx - 1){
      difLimX = dp_var - (xmax - particula.px);
      if (difLimX > min_value){
        nuevaAx = ax - (ps * difLimX - dv * particula.vx);
      }
    }
    return nuevaAx;
}

double colisionLimiteEjeY(Particula& particula, double ymin, int cy, double ymax, double tiempo, double dp_var, int ny, double ay, double ps, double dv) {
    double nuevaY = particula.py + particula.hvy * tiempo;
    double const min_value = 0.0000000001;
    double nuevaAy = ay;
    if (cy == 0){
      double difLimY = dp_var - (particula.py - ymin);
      if (difLimY > min_value){
        nuevaAy = ay + (ps * difLimY - dv * particula.vy);
      }
    }else if(cy == ny - 1){
      double difLimY = dp_var - (ymax - particula.py);
      if (difLimY > min_value){
        nuevaAy = ay - (ps * difLimY - dv * particula.vy);
      }
    }
    return nuevaAy;
}

double colisionLimiteEjeZ(Particula& particula, double zmin, int cz, double zmax, double tiempo, double dp_var, int nz, double az, double ps, double dv) {
    double nuevaZ = particula.pz + particula.hvz * tiempo;
    double const min_value = 0.0000000001;
    double nuevaAz = az;
    if (cz == 0){
      double difLimZ = dp_var - (particula.pz - zmin);
      if (difLimZ > min_value){
        nuevaAz = az + (ps * difLimZ - dv * particula.vz);
      }
    }else if(cz == nz - 1){
      double difLimZ = dp_var - (zmax - particula.pz);
      if (difLimZ > min_value){
        nuevaAz = az - (ps * difLimZ - dv * particula.vz);
      }
    }
    return nuevaAz;
}

void actualizacionColisiones(Particula& particula, double tiempo) {
    double ax = colisionLimiteEjeX()//tengo q ver que le paso
    particula.px = particula.px + particula.hvx * tiempo + ax * (tiempo * tiempo);
    particula.vx = particula.hvx + (ax * tiempo)/2;
    particula.hvx = hvx + ax * tiempo;

    double ay = colisionLimiteEjeY()//tengo q ver que le paso
    particula.py = particula.py + particula.hvy * tiempo + ay * (tiempo * tiempo);
    particula.vy = particula.hvy + (ay * tiempo)/2;
    particula.hvy = hvy + ay * tiempo;

    double az = colisionLimiteEjeZ()//tengo q ver que le paso
    particula.pz = particula.pz + particula.hvz * tiempo + az * (tiempo * tiempo);
    particula.vz = particula.hvz + (az * tiempo)/2;
    particula.hvz = hvz + az * tiempo;

}*/

