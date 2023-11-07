//
// Created by david on 11/7/23.
//

#ifndef ARCOS_BLOCK_H
#define ARCOS_BLOCK_H


class block {
public:
    // Constructor
    block(int sx, int sy, int sz) : sx(sx), sy(sy), sz(sz) { }

    // Destructor
    ~block() {
        // Liberación de recursos o acciones de limpieza en el destructor
    }

private:
    /*Tamaño de bloques*/
    int sx;
    int sy;
    int sz;
};


#endif//ARCOS_BLOCK_H
