#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>

typedef struct
{
    double X, Y, Z, R, Vx, Vy, Vz, Ax, Ay, Az, WQ;

}ParticleOne;

typedef struct
{
    int npar;
    int size;
    int npar_min;
    int npar_max;
    double k_lb, k_dec, k_inc;

    ParticleOne *po;
    ParticleOne *ptr;

} ParticleBundle;

void c_realloc_bundle_(ParticleBundle *pb, int *n)
{
    pb->po = realloc(pb->po, *n * sizeof(ParticleOne));
}

void c_deallocate_bundle_(ParticleBundle *pb)
{
    if (NULL != pb->po)
    {
        free(pb->po);
        pb->po = NULL;
    }
}
