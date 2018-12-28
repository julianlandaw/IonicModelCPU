//
//  Cable1D.h
//
//  Structure for analyzing the dynamics of variables (e.g. nai, ki)
//  as a function of PCL and APD. This structure allows the user to give
//  two APDs, one APD for set number of beats (BEATS1) and the other APD
//  for another set number of beats (BEATS2)
//
//
//  Created by Julian Landaw on 1/7/17.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#ifndef Cable1D_h
#define Cable1D_h

template <template<int> class typecell, int ncells>
class Cable1D
{
public:
    typecell<ncells> Cells;
    double pcl;
    int numstimulus;
    double dx[ncells];
    
    Cable1D();

    Cable1D(double _pcl, int _numstimulus);
    
    void iterate(const int id, double commondt, double t);

    void diffuse(const int id);
    
    void iterateall(double commondt, double t);
    
    void diffuseall();
};

#endif //Cable1D_h
