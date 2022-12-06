//Typically compiled as unit for a larger project
//Alternatively, use `-DMC-SELF-TEST` to create stand alone executable:
// gcc -DMC_SELF_TEST MarchingCubes.c -o mctest; ./mctest

//------------------------------------------------
// MarchingCubes
//------------------------------------------------
//
// MarchingCubes Algorithm
// Version 0.2 - 12/08/2002
//
// Thomas Lewiner thomas.lewiner@polytechnique.org
// Math Dept, PUC-Rio
//
//
// Translated to C by Ziad S. Saad November 30/04
//________________________________________________

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>
#include "LookUpTable.h"
#include "MarchingCubes.h"
#ifdef MC_SELF_TEST
 #include <unistd.h>
#endif
// step size of the arrays of vertices and triangles
#define ALLOC_SIZE 65536

static int debug;
void set_suma_debug(int dbg)
{
   debug = dbg;
   return;
}
//_____________________________________________________________________________
// print cube for debug
void print_cube(MCB *mcb) { printf( "\t%f %f %f %f %f %f %f %f\n", mcb->cube[0], mcb->cube[1], mcb->cube[2], mcb->cube[3], mcb->cube[4], mcb->cube[5], mcb->cube[6], mcb->cube[7]) ; }
//_____________________________________________________________________________

void set_resolution( MCB *mcb, int size_x,  int size_y,  int size_z ) 
{ 
   mcb->size_x = size_x ;  mcb->size_y = size_y ;  mcb->size_z = size_z ; 
   return;
}
void set_method    ( MCB *mcb, int originalMC ) {
    /* originalMC = false is the default */ 
    mcb->originalMC = originalMC ; 
    return;
}

  // Data access
float get_data  (  MCB *mcb, int i,  int j,  int k )  { 
   return (mcb->data[ i + j*mcb->size_x + k*mcb->size_x*mcb->size_y]) ; 
}

void set_data  (  MCB *mcb, float val,  int i,  int j,  int k ) {
  (mcb->data[ i + j*mcb->size_x + k*mcb->size_x*mcb->size_y] = val) ; 
}
int   get_x_vert(  MCB *mcb , int i,  int j,  int k )  { return (mcb->x_verts[ i + j*mcb->size_x + k*mcb->size_x*mcb->size_y] ); }
int   get_y_vert(  MCB *mcb , int i,  int j,  int k )  { return (mcb->y_verts[ i + j*mcb->size_x + k*mcb->size_x*mcb->size_y] ); }
int   get_z_vert(  MCB *mcb , int i,  int j,  int k )  { return (mcb->z_verts[ i + j*mcb->size_x + k*mcb->size_x*mcb->size_y] ); }
void  set_x_vert(  MCB *mcb , int val,  int i,  int j,  int k ) { (mcb->x_verts[ i + j*mcb->size_x + k*mcb->size_x*mcb->size_y] = val ); }
void  set_y_vert(  MCB *mcb , int val,  int i,  int j,  int k ) { (mcb->y_verts[ i + j*mcb->size_x + k*mcb->size_x*mcb->size_y] = val ); }
void  set_z_vert(  MCB *mcb , int val,  int i,  int j,  int k ) { (mcb->z_verts[ i + j*mcb->size_x + k*mcb->size_x*mcb->size_y] = val ); }


//_____________________________________________________________________________
// Constructor
MCB *MarchingCubes( int size_x, int size_y , int size_z  )
{
// defaults are -1 for all size_ -----------------------------------------------------------------------------
  MCB *mcb=NULL;
  mcb = (MCB *)malloc(sizeof(MCB));
  mcb->originalMC = false;
  mcb->size_x = size_x;
  mcb->size_y = size_y;
  mcb->size_z = size_z;
  mcb->data    =  (float*)NULL;
  mcb->x_verts = ( int *)NULL;
  mcb->y_verts = ( int *)NULL;
  mcb->z_verts =  ( int *)NULL;
  mcb->nverts  =  0;
  mcb->ntrigs  =  0;
  mcb->Nverts  =  0;
  mcb->Ntrigs  =  0;
  mcb->vertices = ( Vertex *)NULL;
  mcb->triangles =(Triangle*)NULL;
  mcb->_case = 0;             /* Was uninitialized and causing weird crashes on linux! ZSS: Oct 06 */
  return(mcb);
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Destructor
void FreeMarchingCubes(MCB *mcb)
//-----------------------------------------------------------------------------
{
  clean_all(mcb) ;
  return;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// main algorithm
void run(MCB *mcb)
//-----------------------------------------------------------------------------
{
   int p;
   if (debug) printf("Marching Cubes begin: cpu %f\n", (double)clock() ) ;

   compute_intersection_points( mcb) ;

   for( mcb->k = 0 ; mcb->k < mcb->size_z-1 ; mcb->k++ )
   for( mcb->j = 0 ; mcb->j < mcb->size_y-1 ; mcb->j++ )
   for( mcb->i = 0 ; mcb->i < mcb->size_x-1 ; mcb->i++ )
   {
    mcb->lut_entry = 0 ;
    for(  p = 0 ; p < 8 ; ++p )
    {
      mcb->cube[p] = get_data( mcb, mcb->i+((p^(p>>1))&1), mcb->j+((p>>1)&1), mcb->k+((p>>2)&1) ) ;
      if( fabs( mcb->cube[p] ) < FLT_EPSILON ) mcb->cube[p] = FLT_EPSILON ;
      if( mcb->cube[p] > 0 ) mcb->lut_entry += 1 << p ;
    }
   /*
    if( ( mcb->cube[0] = get_data( mcb, mcb->i , mcb->j , mcb->k ) ) > 0 ) mcb->lut_entry +=   1 ;
    if( ( mcb->cube[1] = get_data(mcb, mcb->i+1, mcb->j , mcb->k ) ) > 0 ) mcb->lut_entry +=   2 ;
    if( ( mcb->cube[2] = get_data(mcb, mcb->i+1,mcb->j+1, mcb->k ) ) > 0 ) mcb->lut_entry +=   4 ;
    if( ( mcb->cube[3] = get_data(mcb,  mcb->i ,mcb->j+1, mcb->k ) ) > 0 ) mcb->lut_entry +=   8 ;
    if( ( mcb->cube[4] = get_data(mcb,  mcb->i , mcb->j ,mcb->k+1) ) > 0 ) mcb->lut_entry +=  16 ;
    if( ( mcb->cube[5] = get_data(mcb, mcb->i+1, mcb->j ,mcb->k+1) ) > 0 ) mcb->lut_entry +=  32 ;
    if( ( mcb->cube[6] = get_data(mcb, mcb->i+1,mcb->j+1,mcb->k+1) ) > 0 ) mcb->lut_entry +=  64 ;
    if( ( mcb->cube[7] = get_data(mcb,  mcb->i ,mcb->j+1,mcb->k+1) ) > 0 ) mcb->lut_entry += 128 ;
   */
    process_cube( mcb) ;
   }

   if (debug) { 
      printf("Marching Cubes end: cpu %g\n", (double)clock() ) ;
      for( mcb->i = 0 ; mcb->i < 15 ; mcb->i++ )
      {
       printf("  %7d cases %d\n", mcb->N[mcb->i], mcb->i ) ;
      }
   }
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// init temporary structures (must set sizes before call)
void init_temps(MCB *mcb)
//-----------------------------------------------------------------------------
{
  mcb->data = (float*)calloc(mcb->size_x * mcb->size_y * mcb->size_z, sizeof(float)) ;
  mcb->x_verts = (int*)calloc(mcb->size_x * mcb->size_y * mcb->size_z, sizeof(int)) ;
  mcb->y_verts = (int*)calloc(mcb->size_x * mcb->size_y * mcb->size_z, sizeof(int));
  mcb->z_verts = (int*)calloc(mcb->size_x * mcb->size_y * mcb->size_z, sizeof(int)) ;

  memset( mcb->x_verts, -1, mcb->size_x * mcb->size_y * mcb->size_z * sizeof( int ) ) ;
  memset( mcb->y_verts, -1, mcb->size_x * mcb->size_y * mcb->size_z * sizeof( int ) ) ;
  memset( mcb->z_verts, -1, mcb->size_x * mcb->size_y * mcb->size_z * sizeof( int ) ) ;

  memset( mcb->N, 0, 15 * sizeof(int) ) ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// init all structures (must set sizes before call)
void init_all (MCB *mcb)
//-----------------------------------------------------------------------------
{
  init_temps(mcb) ;

  mcb->nverts = mcb->ntrigs = 0 ;
  mcb->Nverts = mcb->Ntrigs = ALLOC_SIZE ;
  mcb->vertices  = (Vertex*)calloc(mcb->Nverts, sizeof(Vertex)) ;
  mcb->triangles = (Triangle*)calloc(mcb->Ntrigs, sizeof(Triangle));
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// clean temporary structures
void clean_temps(MCB *mcb)
//-----------------------------------------------------------------------------
{
  free(mcb->data); 
  free(mcb->x_verts);
  free(mcb->y_verts);
  free(mcb->z_verts);

  mcb->data     = (float*)NULL ;
  mcb->x_verts  = (int*)NULL ;
  mcb->y_verts  = (int*)NULL ;
  mcb->z_verts  = (int*)NULL ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// clean all structures
void clean_all(MCB *mcb)
//-----------------------------------------------------------------------------
{
  clean_temps(mcb) ;
  free(mcb->vertices)  ;
  free(mcb->triangles) ;
  mcb->vertices  = (Vertex   *)NULL ;
  mcb->triangles = (Triangle *)NULL ;
  mcb->nverts = mcb->ntrigs = 0 ;
  mcb->Nverts = mcb->Ntrigs = 0 ;

  mcb->size_x = mcb->size_y = mcb->size_z = -1 ;
}
//_____________________________________________________________________________


//_____________________________________________________________________________
// Compute the intersection points
void compute_intersection_points(MCB *mcb )
//-----------------------------------------------------------------------------
{
  for( mcb->k = 0 ; mcb->k < mcb->size_z ; mcb->k++ )
  for( mcb->j = 0 ; mcb->j < mcb->size_y ; mcb->j++ )
  for( mcb->i = 0 ; mcb->i < mcb->size_x ; mcb->i++ )
  {
    mcb->cube[0] = get_data( mcb, mcb->i, mcb->j, mcb->k ) ;
    if( mcb->i < mcb->size_x - 1 ) mcb->cube[1] = get_data(mcb, mcb->i+1, mcb->j , mcb->k ) ;
    else                 mcb->cube[1] = mcb->cube[0] ;

    if( mcb->j < mcb->size_y - 1 ) mcb->cube[3] = get_data( mcb, mcb->i ,mcb->j+1, mcb->k ) ;
    else                 mcb->cube[3] = mcb->cube[0] ;

    if( mcb->k < mcb->size_z - 1 ) mcb->cube[4] = get_data( mcb, mcb->i , mcb->j ,mcb->k+1) ;
    else                 mcb->cube[4] = mcb->cube[0] ;

    if( fabs( mcb->cube[0] ) < FLT_EPSILON ) mcb->cube[0] = FLT_EPSILON ;
    if( fabs( mcb->cube[1] ) < FLT_EPSILON ) mcb->cube[1] = FLT_EPSILON ;
    if( fabs( mcb->cube[3] ) < FLT_EPSILON ) mcb->cube[3] = FLT_EPSILON ;
    if( fabs( mcb->cube[4] ) < FLT_EPSILON ) mcb->cube[4] = FLT_EPSILON ;

    if( mcb->cube[0] < 0 )
    {
      if( mcb->cube[1] > 0 ) set_x_vert( mcb, add_x_vertex( mcb), mcb->i,mcb->j,mcb->k ) ;
      if( mcb->cube[3] > 0 ) set_y_vert( mcb, add_y_vertex( mcb), mcb->i,mcb->j,mcb->k ) ;
      if( mcb->cube[4] > 0 ) set_z_vert( mcb, add_z_vertex( mcb), mcb->i,mcb->j,mcb->k ) ;
    }
    else
    {
      if( mcb->cube[1] < 0 ) set_x_vert( mcb, add_x_vertex( mcb), mcb->i,mcb->j,mcb->k ) ;
      if( mcb->cube[3] < 0 ) set_y_vert( mcb, add_y_vertex( mcb), mcb->i,mcb->j,mcb->k ) ;
      if( mcb->cube[4] < 0 ) set_z_vert( mcb, add_z_vertex( mcb), mcb->i,mcb->j,mcb->k ) ;
    }
  }
}
//_____________________________________________________________________________

//_____________________________________________________________________________
// Test a face
// if face>0 return true if the face contains a part of the surface
int test_face( MCB *mcb, schar face )
//-----------------------------------------------------------------------------
{
  float A,B,C,D ;

  switch( face )
  {
  case -1 : case 1 :  A = mcb->cube[0] ;  B = mcb->cube[4] ;  C = mcb->cube[5] ;  D = mcb->cube[1] ;  break ;
  case -2 : case 2 :  A = mcb->cube[1] ;  B = mcb->cube[5] ;  C = mcb->cube[6] ;  D = mcb->cube[2] ;  break ;
  case -3 : case 3 :  A = mcb->cube[2] ;  B = mcb->cube[6] ;  C = mcb->cube[7] ;  D = mcb->cube[3] ;  break ;
  case -4 : case 4 :  A = mcb->cube[3] ;  B = mcb->cube[7] ;  C = mcb->cube[4] ;  D = mcb->cube[0] ;  break ;
  case -5 : case 5 :  A = mcb->cube[0] ;  B = mcb->cube[3] ;  C = mcb->cube[2] ;  D = mcb->cube[1] ;  break ;
  case -6 : case 6 :  A = mcb->cube[4] ;  B = mcb->cube[7] ;  C = mcb->cube[6] ;  D = mcb->cube[5] ;  break ;
  default : printf( "Invalid face code %d %d %d: %d\n",  mcb->i,  mcb->j,  mcb->k, face ) ;  print_cube(mcb) ;  A = B = C = D = 0 ;
  };

  if( fabs( A*C - B*D ) < FLT_EPSILON )
    return face >= 0 ;
  return face * A * ( A*C - B*D ) >= 0  ;  // face and A invert signs
}

//_____________________________________________________________________________
// Test the interior of a cube
// if s == 7, return true  if the interior is empty
// if s ==-7, return false if the interior is empty
int test_interior( MCB *mcb, schar s )
//-----------------------------------------------------------------------------
{
  float t, At=0, Bt=0, Ct=0, Dt=0, a, b ;
  char  test =  0 ;
  char  edge = -1 ; // reference edge of the triangulation

  switch( mcb->_case )
  {
  case  4 :
  case 10 :
    a = ( mcb->cube[4] - mcb->cube[0] ) * ( mcb->cube[6] - mcb->cube[2] ) - ( mcb->cube[7] - mcb->cube[3] ) * ( mcb->cube[5] - mcb->cube[1] ) ;
    b =  mcb->cube[2] * ( mcb->cube[4] - mcb->cube[0] ) + mcb->cube[0] * ( mcb->cube[6] - mcb->cube[2] )
             - mcb->cube[1] * ( mcb->cube[7] - mcb->cube[3] ) - mcb->cube[3] * ( mcb->cube[5] - mcb->cube[1] ) ;
    t = - b / (2*a) ;
    if( t<0 || t>1 ) return s>0 ;

    At = mcb->cube[0] + ( mcb->cube[4] - mcb->cube[0] ) * t ;
    Bt = mcb->cube[3] + ( mcb->cube[7] - mcb->cube[3] ) * t ;
    Ct = mcb->cube[2] + ( mcb->cube[6] - mcb->cube[2] ) * t ;
    Dt = mcb->cube[1] + ( mcb->cube[5] - mcb->cube[1] ) * t ;
    break ;

  case  6 :
  case  7 :
  case 12 :
  case 13 :
    switch( mcb->_case )
    {
    case  6 : edge = test6 [mcb->config][2] ; break ;
    case  7 : edge = test7 [mcb->config][4] ; break ;
    case 12 : edge = test12[mcb->config][3] ; break ;
    case 13 : edge = tiling13_5_1[mcb->config][mcb->subconfig][0] ; break ;
    }
    switch( edge )
    {
    case  0 :
      t  = mcb->cube[0] / ( mcb->cube[0] - mcb->cube[1] ) ;
      At = 0 ;
      Bt = mcb->cube[3] + ( mcb->cube[2] - mcb->cube[3] ) * t ;
      Ct = mcb->cube[7] + ( mcb->cube[6] - mcb->cube[7] ) * t ;
      Dt = mcb->cube[4] + ( mcb->cube[5] - mcb->cube[4] ) * t ;
      break ;
    case  1 :
      t  = mcb->cube[1] / ( mcb->cube[1] - mcb->cube[2] ) ;
      At = 0 ;
      Bt = mcb->cube[0] + ( mcb->cube[3] - mcb->cube[0] ) * t ;
      Ct = mcb->cube[4] + ( mcb->cube[7] - mcb->cube[4] ) * t ;
      Dt = mcb->cube[5] + ( mcb->cube[6] - mcb->cube[5] ) * t ;
      break ;
    case  2 :
      t  = mcb->cube[2] / ( mcb->cube[2] - mcb->cube[3] ) ;
      At = 0 ;
      Bt = mcb->cube[1] + ( mcb->cube[0] - mcb->cube[1] ) * t ;
      Ct = mcb->cube[5] + ( mcb->cube[4] - mcb->cube[5] ) * t ;
      Dt = mcb->cube[6] + ( mcb->cube[7] - mcb->cube[6] ) * t ;
      break ;
    case  3 :
      t  = mcb->cube[3] / ( mcb->cube[3] - mcb->cube[0] ) ;
      At = 0 ;
      Bt = mcb->cube[2] + ( mcb->cube[1] - mcb->cube[2] ) * t ;
      Ct = mcb->cube[6] + ( mcb->cube[5] - mcb->cube[6] ) * t ;
      Dt = mcb->cube[7] + ( mcb->cube[4] - mcb->cube[7] ) * t ;
      break ;
    case  4 :
      t  = mcb->cube[4] / ( mcb->cube[4] - mcb->cube[5] ) ;
      At = 0 ;
      Bt = mcb->cube[7] + ( mcb->cube[6] - mcb->cube[7] ) * t ;
      Ct = mcb->cube[3] + ( mcb->cube[2] - mcb->cube[3] ) * t ;
      Dt = mcb->cube[0] + ( mcb->cube[1] - mcb->cube[0] ) * t ;
      break ;
    case  5 :
      t  = mcb->cube[5] / ( mcb->cube[5] - mcb->cube[6] ) ;
      At = 0 ;
      Bt = mcb->cube[4] + ( mcb->cube[7] - mcb->cube[4] ) * t ;
      Ct = mcb->cube[0] + ( mcb->cube[3] - mcb->cube[0] ) * t ;
      Dt = mcb->cube[1] + ( mcb->cube[2] - mcb->cube[1] ) * t ;
      break ;
    case  6 :
      t  = mcb->cube[6] / ( mcb->cube[6] - mcb->cube[7] ) ;
      At = 0 ;
      Bt = mcb->cube[5] + ( mcb->cube[4] - mcb->cube[5] ) * t ;
      Ct = mcb->cube[1] + ( mcb->cube[0] - mcb->cube[1] ) * t ;
      Dt = mcb->cube[2] + ( mcb->cube[3] - mcb->cube[2] ) * t ;
      break ;
    case  7 :
      t  = mcb->cube[7] / ( mcb->cube[7] - mcb->cube[4] ) ;
      At = 0 ;
      Bt = mcb->cube[6] + ( mcb->cube[5] - mcb->cube[6] ) * t ;
      Ct = mcb->cube[2] + ( mcb->cube[1] - mcb->cube[2] ) * t ;
      Dt = mcb->cube[3] + ( mcb->cube[0] - mcb->cube[3] ) * t ;
      break ;
    case  8 :
      t  = mcb->cube[0] / ( mcb->cube[0] - mcb->cube[4] ) ;
      At = 0 ;
      Bt = mcb->cube[3] + ( mcb->cube[7] - mcb->cube[3] ) * t ;
      Ct = mcb->cube[2] + ( mcb->cube[6] - mcb->cube[2] ) * t ;
      Dt = mcb->cube[1] + ( mcb->cube[5] - mcb->cube[1] ) * t ;
      break ;
    case  9 :
      t  = mcb->cube[1] / ( mcb->cube[1] - mcb->cube[5] ) ;
      At = 0 ;
      Bt = mcb->cube[0] + ( mcb->cube[4] - mcb->cube[0] ) * t ;
      Ct = mcb->cube[3] + ( mcb->cube[7] - mcb->cube[3] ) * t ;
      Dt = mcb->cube[2] + ( mcb->cube[6] - mcb->cube[2] ) * t ;
      break ;
    case 10 :
      t  = mcb->cube[2] / ( mcb->cube[2] - mcb->cube[6] ) ;
      At = 0 ;
      Bt = mcb->cube[1] + ( mcb->cube[5] - mcb->cube[1] ) * t ;
      Ct = mcb->cube[0] + ( mcb->cube[4] - mcb->cube[0] ) * t ;
      Dt = mcb->cube[3] + ( mcb->cube[7] - mcb->cube[3] ) * t ;
      break ;
    case 11 :
      t  = mcb->cube[3] / ( mcb->cube[3] - mcb->cube[7] ) ;
      At = 0 ;
      Bt = mcb->cube[2] + ( mcb->cube[6] - mcb->cube[2] ) * t ;
      Ct = mcb->cube[1] + ( mcb->cube[5] - mcb->cube[1] ) * t ;
      Dt = mcb->cube[0] + ( mcb->cube[4] - mcb->cube[0] ) * t ;
      break ;
    default : printf( "Invalid edge %d\n", edge ) ;  print_cube(mcb) ;  break ;
    }
    break ;

  default : printf( "Invalid ambiguous case %d\n", mcb->_case ) ;  print_cube(mcb) ;  break ;
  }

  if( At >= 0 ) test ++ ;
  if( Bt >= 0 ) test += 2 ;
  if( Ct >= 0 ) test += 4 ;
  if( Dt >= 0 ) test += 8 ;
  switch( test )
  {
  case  0 : return s>0 ;
  case  1 : return s>0 ;
  case  2 : return s>0 ;
  case  3 : return s>0 ;
  case  4 : return s>0 ;
  case  5 : if( At * Ct - Bt * Dt <  FLT_EPSILON ) return s>0 ; break ;
  case  6 : return s>0 ;
  case  7 : return s<0 ;
  case  8 : return s>0 ;
  case  9 : return s>0 ;
  case 10 : if( At * Ct - Bt * Dt >= FLT_EPSILON ) return s>0 ; break ;
  case 11 : return s<0 ;
  case 12 : return s>0 ;
  case 13 : return s<0 ;
  case 14 : return s<0 ;
  case 15 : return s<0 ;
  }

  return s<0 ;
}
//_____________________________________________________________________________

//_____________________________________________________________________________
// Process a unit cube
void process_cube( MCB *mcb)
//-----------------------------------------------------------------------------
{
  int   v12 = -1 ;
  /* print_cube(mcb) ; 
  fprintf (stderr,"_case=%d\n", mcb->_case);
  fprintf (stderr,"N=%d\n", mcb->N[mcb->_case]);*/
  if (mcb->_case >= N_MAX) {
   fprintf (stderr,"Unexpected _case value of %d\nResetting to 0.\n",mcb->_case);
   mcb->_case = 0; 
  }
  mcb->N[mcb->_case]++ ;

  if( mcb->originalMC )
  {
    char nt = 0 ;
    while( casesClassic[mcb->lut_entry][3*nt] != -1 ) nt++ ;
    add_triangle(mcb, casesClassic[mcb->lut_entry], nt, -1 ) ;
    return ;
  }

  mcb->_case   = cases[mcb->lut_entry][0] ;
  mcb->config = cases[mcb->lut_entry][1] ;
  mcb->subconfig = 0 ;

  switch( mcb->_case )
  {
  case  0 :
    break ;

  case  1 :
    add_triangle(mcb, tiling1[mcb->config], 1, -1) ;
    break ;

  case  2 :
    add_triangle(mcb, tiling2[mcb->config], 2, -1) ;
    break ;

  case  3 :
    if( test_face(mcb, test3[mcb->config]) )
      add_triangle(mcb,  tiling3_2[mcb->config], 4, -1) ; // 3.2
    else
      add_triangle(mcb,  tiling3_1[mcb->config], 2, -1) ; // 3.1
    break ;

  case  4 :
    if( test_interior(mcb, test4[mcb->config]) )
      add_triangle(mcb, tiling4_1[mcb->config], 2, -1) ; // 4.1.1
    else
      add_triangle(mcb, tiling4_2[mcb->config], 6, -1) ; // 4.1.2
    break ;

  case  5 :
    add_triangle(mcb, tiling5[mcb->config], 3, -1) ;
    break ;

  case  6 :
    if( test_face(mcb, test6[mcb->config][0]) )
      add_triangle(mcb, tiling6_2[mcb->config], 5, -1) ; // 6.2
    else
    {
      if( test_interior(mcb, test6[mcb->config][1]) )
        add_triangle(mcb, tiling6_1_1[mcb->config], 3, -1) ; // 6.1.1
      else
    {
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling6_1_2[mcb->config], 9 , v12) ; // 6.1.2
      }
    }
    break ;
  case  7 :
    if( test_face(mcb, test7[mcb->config][0] ) ) mcb->subconfig +=  1 ;
    if( test_face(mcb, test7[mcb->config][1] ) ) mcb->subconfig +=  2 ;
    if( test_face(mcb, test7[mcb->config][2] ) ) mcb->subconfig +=  4 ;
    switch( mcb->subconfig )
      {
      case 0 :
        add_triangle(mcb, tiling7_1[mcb->config], 3, -1) ; break ;
      case 1 :
        add_triangle(mcb, tiling7_2[mcb->config][0], 5, -1) ; break ;
      case 2 :
        add_triangle(mcb, tiling7_2[mcb->config][1], 5, -1) ; break ;
      case 3 :
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling7_3[mcb->config][0], 9, v12 ) ; break ;
      case 4 :
        add_triangle(mcb, tiling7_2[mcb->config][2], 5, -1) ; break ;
      case 5 :
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling7_3[mcb->config][1], 9, v12 ) ; break ;
      case 6 :
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling7_3[mcb->config][2], 9, v12 ) ; break ;
      case 7 :
        if( test_interior(mcb, test7[mcb->config][3]) )
          add_triangle(mcb, tiling7_4_2[mcb->config], 9, -1) ;
        else
          add_triangle(mcb, tiling7_4_1[mcb->config], 5, -1) ;
        break ;
      };
    break ;

  case  8 :
    add_triangle(mcb, tiling8[mcb->config], 2, -1) ;
    break ;

  case  9 :
    add_triangle(mcb, tiling9[mcb->config], 4, -1) ;
    break ;

  case 10 :
    if( test_face(mcb, test10[mcb->config][0]) )
    {
      if( test_face(mcb, test10[mcb->config][1]) )
        add_triangle(mcb, tiling10_1_1_[mcb->config], 4, -1) ; // 10.1.1
      else
      {
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling10_2[mcb->config], 8, v12 ) ; // 10.2
      }
    }
    else
    {
      if( test_face(mcb, test10[mcb->config][1]) )
      {
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling10_2_[mcb->config], 8, v12 ) ; // 10.2
      }
      else
      {
        if( test_interior(mcb, test10[mcb->config][2]) )
          add_triangle(mcb, tiling10_1_1[mcb->config], 4, -1) ; // 10.1.1
        else
          add_triangle(mcb, tiling10_1_2[mcb->config], 8, -1) ; // 10.1.2
      }
    }
    break ;

  case 11 :
    add_triangle(mcb, tiling11[mcb->config], 4, -1) ;
    break ;

  case 12 :
    if( test_face(mcb, test12[mcb->config][0]) )
    {
      if( test_face(mcb, test12[mcb->config][1]) )
        add_triangle(mcb, tiling12_1_1_[mcb->config], 4, -1) ; // 12.1.1
      else
      {
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling12_2[mcb->config], 8, v12 ) ; // 12.2
      }
    }
    else
    {
      if( test_face(mcb, test12[mcb->config][1]) )
      {
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling12_2_[mcb->config], 8, v12 ) ; // 12.2
      }
      else
      {
        if( test_interior(mcb, test12[mcb->config][2]) )
          add_triangle(mcb, tiling12_1_1[mcb->config], 4, -1) ; // 12.1.1
        else
          add_triangle(mcb, tiling12_1_2[mcb->config], 8, -1) ; // 12.1.2
      }
    }
    break ;

  case 13 :
    if( test_face(mcb,  test13[mcb->config][0] ) ) mcb->subconfig +=  1 ;
    if( test_face(mcb,  test13[mcb->config][1] ) ) mcb->subconfig +=  2 ;
    if( test_face(mcb,  test13[mcb->config][2] ) ) mcb->subconfig +=  4 ;
    if( test_face(mcb,  test13[mcb->config][3] ) ) mcb->subconfig +=  8 ;
    if( test_face(mcb,  test13[mcb->config][4] ) ) mcb->subconfig += 16 ;
    if( test_face(mcb,  test13[mcb->config][5] ) ) mcb->subconfig += 32 ;
    switch( subconfig13[mcb->subconfig] )
    {
      case 0 :/* 13.1 */
        add_triangle(mcb,  tiling13_1[mcb->config], 4, -1) ; break ;

      case 1 :/* 13.2 */
        add_triangle(mcb,  tiling13_2[mcb->config][0], 6, -1) ; break ;
      case 2 :/* 13.2 */
        add_triangle(mcb,  tiling13_2[mcb->config][1], 6, -1) ; break ;
      case 3 :/* 13.2 */
        add_triangle(mcb,  tiling13_2[mcb->config][2], 6, -1) ; break ;
      case 4 :/* 13.2 */
        add_triangle(mcb,  tiling13_2[mcb->config][3], 6, -1) ; break ;
      case 5 :/* 13.2 */
        add_triangle(mcb,  tiling13_2[mcb->config][4], 6, -1) ; break ;
      case 6 :/* 13.2 */
        add_triangle(mcb,  tiling13_2[mcb->config][5], 6, -1) ; break ;

      case 7 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb,  tiling13_3[mcb->config][0], 10, v12 ) ; break ;
      case 8 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb,  tiling13_3[mcb->config][1], 10, v12 ) ; break ;
      case 9 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb,  tiling13_3[mcb->config][2], 10, v12 ) ; break ;
      case 10 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb,  tiling13_3[mcb->config][3], 10, v12 ) ; break ;
      case 11 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb,  tiling13_3[mcb->config][4], 10, v12 ) ; break ;
      case 12 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb,  tiling13_3[mcb->config][5], 10, v12 ) ; break ;
      case 13 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb,  tiling13_3[mcb->config][6], 10, v12 ) ; break ;
      case 14 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb,  tiling13_3[mcb->config][7], 10, v12 ) ; break ;
      case 15 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb,  tiling13_3[mcb->config][8], 10, v12 ) ; break ;
      case 16 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb,  tiling13_3[mcb->config][9], 10, v12 ) ; break ;
      case 17 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb,  tiling13_3[mcb->config][10], 10, v12 ) ; break ;
      case 18 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb,  tiling13_3[mcb->config][11], 10, v12 ) ; break ;

      case 19 :/* 13.4 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_4[mcb->config][0], 12, v12 ) ; break ;
      case 20 :/* 13.4 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_4[mcb->config][1], 12, v12 ) ; break ;
      case 21 :/* 13.4 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_4[mcb->config][2], 12, v12 ) ; break ;
      case 22 :/* 13.4 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_4[mcb->config][3], 12, v12 ) ; break ;

      case 23 :/* 13.5 */
        mcb->subconfig = 0 ;
        if( test_interior(mcb, test13[mcb->config][6] ) )
          add_triangle(mcb, tiling13_5_1[mcb->config][0], 6, -1) ;
        else
          add_triangle(mcb, tiling13_5_2[mcb->config][0], 10, -1) ;
        break ;
      case 24 :/* 13.5 */
        mcb->subconfig = 1 ;
        if( test_interior(mcb, test13[mcb->config][6] ) )
          add_triangle(mcb, tiling13_5_1[mcb->config][1], 6, -1) ;
        else
          add_triangle(mcb, tiling13_5_2[mcb->config][1], 10, -1) ;
        break ;
      case 25 :/* 13.5 */
        mcb->subconfig = 2 ;
        if( test_interior(mcb, test13[mcb->config][6] ) )
          add_triangle(mcb, tiling13_5_1[mcb->config][2], 6, -1) ;
        else
          add_triangle(mcb, tiling13_5_2[mcb->config][2], 10, -1) ;
        break ;
      case 26 :/* 13.5 */
        mcb->subconfig = 3 ;
        if( test_interior(mcb, test13[mcb->config][6] ) )
          add_triangle(mcb, tiling13_5_1[mcb->config][3], 6, -1) ;
        else
          add_triangle(mcb, tiling13_5_2[mcb->config][3], 10, -1) ;
        break ;


      case 27 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_3_[mcb->config][0], 10, v12 ) ; break ;
      case 28 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_3_[mcb->config][1], 10, v12 ) ; break ;
      case 29 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_3_[mcb->config][2], 10, v12 ) ; break ;
      case 30 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_3_[mcb->config][3], 10, v12 ) ; break ;
      case 31 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_3_[mcb->config][4], 10, v12 ) ; break ;
      case 32 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_3_[mcb->config][5], 10, v12 ) ; break ;
      case 33 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_3_[mcb->config][6], 10, v12 ) ; break ;
      case 34 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_3_[mcb->config][7], 10, v12 ) ; break ;
      case 35 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_3_[mcb->config][8], 10, v12 ) ; break ;
      case 36 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_3_[mcb->config][9], 10, v12 ) ; break ;
      case 37 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_3_[mcb->config][10], 10, v12 ) ; break ;
      case 38 :/* 13.3 */
        v12 = add_c_vertex(mcb) ;
        add_triangle(mcb, tiling13_3_[mcb->config][11], 10, v12 ) ; break ;

      case 39 :/* 13.2 */
        add_triangle(mcb, tiling13_2_[mcb->config][0], 6, -1) ; break ;
      case 40 :/* 13.2 */
        add_triangle(mcb, tiling13_2_[mcb->config][1], 6, -1) ; break ;
      case 41 :/* 13.2 */
        add_triangle(mcb, tiling13_2_[mcb->config][2], 6, -1) ; break ;
      case 42 :/* 13.2 */
        add_triangle(mcb, tiling13_2_[mcb->config][3], 6, -1) ; break ;
      case 43 :/* 13.2 */
        add_triangle(mcb, tiling13_2_[mcb->config][4], 6, -1) ; break ;
      case 44 :/* 13.2 */
        add_triangle(mcb, tiling13_2_[mcb->config][5], 6, -1) ; break ;

      case 45 :/* 13.1 */
        add_triangle(mcb, tiling13_1_[mcb->config], 4, -1) ; break ;

      default :
        printf("Marching Cubes: Impossible case 13?\n" ) ;  print_cube(mcb) ;
      }
      break ;

  case 14 :
    add_triangle(mcb, tiling14[mcb->config], 4, -1) ;
    break ;
  };
}

//_____________________________________________________________________________



//_____________________________________________________________________________
// Adding triangles
void add_triangle( MCB *mcb , const char* trig, char n, int v12 )
//-----------------------------------------------------------------------------
{
  int   t, tv[3] ;
//printf( "+>> %d %d %d\n", mcb->i  , mcb->j , mcb->k);

  for( t = 0 ; t < 3*n ; t++ )
  {
    switch( trig[t] )
    {
    case  0 : tv[ t % 3 ] = get_x_vert(mcb,  mcb->i  , mcb->j , mcb->k ) ; break ;
    case  1 : tv[ t % 3 ] = get_y_vert(mcb, mcb->i +1, mcb->j , mcb->k ) ; break ;
    case  2 : tv[ t % 3 ] = get_x_vert(mcb,  mcb->i  ,mcb->j+1, mcb->k ) ; break ;
    case  3 : tv[ t % 3 ] = get_y_vert(mcb,  mcb->i  , mcb->j , mcb->k ) ; break ;
    case  4 : tv[ t % 3 ] = get_x_vert(mcb,  mcb->i  , mcb->j ,mcb->k+1) ; break ;
    case  5 : tv[ t % 3 ] = get_y_vert(mcb, mcb->i +1, mcb->j ,mcb->k+1) ; break ;
    case  6 : tv[ t % 3 ] = get_x_vert(mcb,  mcb->i  ,mcb->j+1,mcb->k+1) ; break ;
    case  7 : tv[ t % 3 ] = get_y_vert(mcb,  mcb->i  , mcb->j ,mcb->k+1) ; break ;
    case  8 : tv[ t % 3 ] = get_z_vert(mcb,  mcb->i  , mcb->j , mcb->k ) ; break ;
    case  9 : tv[ t % 3 ] = get_z_vert(mcb, mcb->i +1, mcb->j , mcb->k ) ; break ;
    case 10 : tv[ t % 3 ] = get_z_vert(mcb, mcb->i +1,mcb->j+1, mcb->k ) ; break ;
    case 11 : tv[ t % 3 ] = get_z_vert(mcb,  mcb->i  ,mcb->j+1, mcb->k ) ; break ;
    case 12 : tv[ t % 3 ] = v12 ; break ;
    default : break ;
    }

    if( tv[t%3] == -1 )
    {
      printf("Marching Cubes: invalid triangle %d\n", mcb->ntrigs+1) ;
      print_cube(mcb) ;
    }

    if( t%3 == 2 )
    { 
      Triangle *T = NULL;
      if( mcb->ntrigs >= mcb->Ntrigs )
      {
        Triangle *temp = mcb->triangles ;
        mcb->triangles = (Triangle*)malloc(2*mcb->Ntrigs * sizeof(Triangle));
        memcpy( mcb->triangles, temp, mcb->Ntrigs*sizeof(Triangle) ) ;
        free(temp) ; temp = NULL;
        if (debug) printf("%d allocated triangles\n", mcb->Ntrigs) ;
        mcb->Ntrigs *= 2 ;
      }

      T = mcb->triangles + mcb->ntrigs++ ;
      T->v1    = tv[0] ;
      T->v2    = tv[1] ;
      T->v3    = tv[2] ;
    }
  }
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Calculating gradient

float get_x_grad( MCB *mcb,  int i,  int j,  int k ) 
//-----------------------------------------------------------------------------
{
  if( i > 0 )
  {
    if ( i < mcb->size_x - 1 )
      return (( get_data( mcb, i+1, j, k ) - get_data( mcb, i-1, j, k ) ) / 2) ;
    else
      return (get_data( mcb, i, j, k ) - get_data( mcb, i-1, j, k )) ;
  }
  else
    return (get_data( mcb, i+1, j, k ) - get_data( mcb, i, j, k )) ;
}
//-----------------------------------------------------------------------------

float get_y_grad( MCB *mcb,  int i,  int j,  int k ) 
//-----------------------------------------------------------------------------
{
  if( j > 0 )
  {
    if ( j < mcb->size_y - 1 )
      return (( get_data( mcb, i, j+1, k ) - get_data( mcb, i, j-1, k ) ) / 2) ;
    else
      return (get_data( mcb, i, j, k ) - get_data( mcb, i, j-1, k )) ;
  }
  else
    return (get_data(mcb,  i, j+1, k ) - get_data(mcb, i, j, k )) ;
}
//-----------------------------------------------------------------------------

float get_z_grad( MCB *mcb, int i,  int j,  int k ) 
//-----------------------------------------------------------------------------
{
  if( k > 0 )
  {
    if ( k < mcb->size_z - 1 )
      return (( get_data( mcb, i, j, k+1 ) - get_data( mcb, i, j, k-1 ) ) / 2) ;
    else
      return (get_data( mcb, i, j, k ) - get_data( mcb, i, j, k-1 )) ;
  }
  else
    return (get_data( mcb, i, j, k+1 ) - get_data( mcb, i, j, k )) ;
}
//_____________________________________________________________________________


//_____________________________________________________________________________
// Adding vertices

void test_vertex_addition(MCB *mcb)
{
  if( mcb->nverts >= mcb->Nverts )
  {
    Vertex *temp = mcb->vertices ;
    mcb->vertices =  (Vertex*)malloc(mcb->Nverts*2 * sizeof(Vertex)) ;
    memcpy( mcb->vertices, temp, mcb->Nverts*sizeof(Vertex) ) ;
    free(temp); temp = NULL;
    if (debug) printf("%d allocated vertices\n", mcb->Nverts) ;
    mcb->Nverts *= 2 ;
  }
}


int add_x_vertex(MCB *mcb )
//-----------------------------------------------------------------------------
{
   Vertex *vert;
   float   u;
  
  test_vertex_addition(mcb) ;
  vert = mcb->vertices + mcb->nverts++ ;

  u = ( mcb->cube[0] ) / ( mcb->cube[0] - mcb->cube[1] ) ;

  vert->x      = (float)mcb->i+u;
  vert->y      = (float) mcb->j ;
  vert->z      = (float) mcb->k ;

  vert->nx = (1-u)*get_x_grad(mcb, mcb->i,mcb->j,mcb->k) + u*get_x_grad(mcb, mcb->i+1,mcb->j,mcb->k) ;
  vert->ny = (1-u)*get_y_grad(mcb, mcb->i,mcb->j,mcb->k) + u*get_y_grad(mcb, mcb->i+1,mcb->j,mcb->k) ;
  vert->nz = (1-u)*get_z_grad(mcb, mcb->i,mcb->j,mcb->k) + u*get_z_grad(mcb, mcb->i+1,mcb->j,mcb->k) ;

  u = (float) sqrt( vert->nx * vert->nx + vert->ny * vert->ny +vert->nz * vert->nz ) ;
  
  if( u > 0 )
  {
    vert->nx /= u ;
    vert->ny /= u ;
    vert->nz /= u ;
  }


  return (mcb->nverts-1) ;
}
//-----------------------------------------------------------------------------

int add_y_vertex( MCB *mcb)
//-----------------------------------------------------------------------------
{  Vertex *vert;
   float u;
  test_vertex_addition(mcb) ;
  vert = mcb->vertices + mcb->nverts++ ;

  u = ( mcb->cube[0] ) / ( mcb->cube[0] - mcb->cube[3] ) ;

  vert->x      = (float) mcb->i ;
  vert->y      = (float)mcb->j+u;
  vert->z      = (float) mcb->k ;

  vert->nx = (1-u)*get_x_grad(mcb, mcb->i,mcb->j,mcb->k) + u*get_x_grad(mcb, mcb->i,mcb->j+1,mcb->k) ;
  vert->ny = (1-u)*get_y_grad(mcb, mcb->i,mcb->j,mcb->k) + u*get_y_grad(mcb, mcb->i,mcb->j+1,mcb->k) ;
  vert->nz = (1-u)*get_z_grad(mcb, mcb->i,mcb->j,mcb->k) + u*get_z_grad(mcb, mcb->i,mcb->j+1,mcb->k) ;

  u = (float) sqrt( vert->nx * vert->nx + vert->ny * vert->ny +vert->nz * vert->nz ) ;
  if( u > 0 )
  {
    vert->nx /= u ;
    vert->ny /= u ;
    vert->nz /= u ;
  }

  return (mcb->nverts-1) ;
}
//-----------------------------------------------------------------------------

int add_z_vertex(MCB *mcb )
//-----------------------------------------------------------------------------
{  Vertex *vert;
   float u;
  test_vertex_addition(mcb) ;
  vert = mcb->vertices + mcb->nverts++ ;

  u = ( mcb->cube[0] ) / ( mcb->cube[0] - mcb->cube[4] ) ;

  vert->x      = (float) mcb->i ;
  vert->y      = (float) mcb->j ;
  vert->z      = (float)mcb->k+u;

  vert->nx = (1-u)*get_x_grad(mcb, mcb->i,mcb->j,mcb->k) + u*get_x_grad(mcb, mcb->i,mcb->j,mcb->k+1) ;
  vert->ny = (1-u)*get_y_grad(mcb, mcb->i,mcb->j,mcb->k) + u*get_y_grad(mcb, mcb->i,mcb->j,mcb->k+1) ;
  vert->nz = (1-u)*get_z_grad(mcb, mcb->i,mcb->j,mcb->k) + u*get_z_grad(mcb, mcb->i,mcb->j,mcb->k+1) ;

  u = (float) sqrt( vert->nx * vert->nx + vert->ny * vert->ny +vert->nz * vert->nz ) ;
  if( u > 0 )
  {
    vert->nx /= u ;
    vert->ny /= u ;
    vert->nz /= u ;
  }

  return (mcb->nverts-1) ;
}


int add_c_vertex( MCB *mcb)
//-----------------------------------------------------------------------------
{  Vertex *vert, v;
   float u;
   int   vid ;
  test_vertex_addition(mcb) ;
  vert = mcb->vertices + mcb->nverts++ ;

  u = 0 ;

  vert->x = vert->y = vert->z =  vert->nx = vert->ny = vert->nz = 0 ;

  // Computes the average of the intersection points of the cube
  vid = get_x_vert( mcb, mcb->i , mcb->j , mcb->k ) ;
  if( vid != -1 ) { ++u ;   v = mcb->vertices[vid] ; vert->x += v.x ;  vert->y += v.y ;  vert->z += v.z ;  vert->nx += v.nx ; vert->ny += v.ny ; vert->nz += v.nz ; }
  vid = get_y_vert(mcb, mcb->i+1, mcb->j , mcb->k ) ;
  if( vid != -1 ) { ++u ;   v = mcb->vertices[vid] ; vert->x += v.x ;  vert->y += v.y ;  vert->z += v.z ;  vert->nx += v.nx ; vert->ny += v.ny ; vert->nz += v.nz ; }
  vid = get_x_vert( mcb, mcb->i ,mcb->j+1, mcb->k ) ;
  if( vid != -1 ) { ++u ;   v = mcb->vertices[vid] ; vert->x += v.x ;  vert->y += v.y ;  vert->z += v.z ;  vert->nx += v.nx ; vert->ny += v.ny ; vert->nz += v.nz ; }
  vid = get_y_vert( mcb, mcb->i , mcb->j , mcb->k ) ;
  if( vid != -1 ) { ++u ;   v = mcb->vertices[vid] ; vert->x += v.x ;  vert->y += v.y ;  vert->z += v.z ;  vert->nx += v.nx ; vert->ny += v.ny ; vert->nz += v.nz ; }
  vid = get_x_vert( mcb, mcb->i , mcb->j ,mcb->k+1) ;
  if( vid != -1 ) { ++u ;   v = mcb->vertices[vid] ; vert->x += v.x ;  vert->y += v.y ;  vert->z += v.z ;  vert->nx += v.nx ; vert->ny += v.ny ; vert->nz += v.nz ; }
  vid = get_y_vert(mcb, mcb->i+1, mcb->j ,mcb->k+1) ;
  if( vid != -1 ) { ++u ;   v = mcb->vertices[vid] ; vert->x += v.x ;  vert->y += v.y ;  vert->z += v.z ;  vert->nx += v.nx ; vert->ny += v.ny ; vert->nz += v.nz ; }
  vid = get_x_vert( mcb, mcb->i ,mcb->j+1,mcb->k+1) ;
  if( vid != -1 ) { ++u ;   v = mcb->vertices[vid] ; vert->x += v.x ;  vert->y += v.y ;  vert->z += v.z ;  vert->nx += v.nx ; vert->ny += v.ny ; vert->nz += v.nz ; }
  vid = get_y_vert(mcb,  mcb->i , mcb->j ,mcb->k+1) ;
  if( vid != -1 ) { ++u ;   v = mcb->vertices[vid] ; vert->x += v.x ;  vert->y += v.y ;  vert->z += v.z ;  vert->nx += v.nx ; vert->ny += v.ny ; vert->nz += v.nz ; }
  vid = get_z_vert( mcb, mcb->i , mcb->j , mcb->k ) ;
  if( vid != -1 ) { ++u ;   v = mcb->vertices[vid] ; vert->x += v.x ;  vert->y += v.y ;  vert->z += v.z ;  vert->nx += v.nx ; vert->ny += v.ny ; vert->nz += v.nz ; }
  vid = get_z_vert(mcb, mcb->i+1, mcb->j , mcb->k ) ;
  if( vid != -1 ) { ++u ;  v = mcb->vertices[vid] ; vert->x += v.x ;  vert->y += v.y ;  vert->z += v.z ;  vert->nx += v.nx ; vert->ny += v.ny ; vert->nz += v.nz ; }
  vid = get_z_vert(mcb, mcb->i+1,mcb->j+1, mcb->k ) ;
  if( vid != -1 ) { ++u ;   v = mcb->vertices[vid] ; vert->x += v.x ;  vert->y += v.y ;  vert->z += v.z ;  vert->nx += v.nx ; vert->ny += v.ny ; vert->nz += v.nz ; }
  vid = get_z_vert( mcb, mcb->i ,mcb->j+1, mcb->k ) ;
  if( vid != -1 ) { ++u ;   v = mcb->vertices[vid] ; vert->x += v.x ;  vert->y += v.y ;  vert->z += v.z ;  vert->nx += v.nx ; vert->ny += v.ny ; vert->nz += v.nz ; }

  vert->x  /= u ;
  vert->y  /= u ;
  vert->z  /= u ;

  u = (float) sqrt( vert->nx * vert->nx + vert->ny * vert->ny +vert->nz * vert->nz ) ;
  if( u > 0 )
  {
    vert->nx /= u ;
    vert->ny /= u ;
    vert->nz /= u ;
  }

  return (mcb->nverts-1) ;
}

#ifdef NII2MESH
int marchingCubes(float * img, size_t dim[3], int lo[3], int hi[3], int originalMC, float isolevel, vec3d **vs, vec3i **ts, int *nv, int *nt) {
  MCB * mcp = MarchingCubes(-1, -1, -1);
  int NX = hi[0] - lo[0] + 1;
  int NY = hi[1] - lo[1] + 1;
  int NZ = hi[2] - lo[2] + 1;
  set_resolution( mcp, NX, NY, NZ) ;
  init_all(mcp) ;
  float * im = mcp->data;
  int i = 0;
  int inX = dim[0];
  int inXY = dim[0] * dim[1];
  for (int z=0;z<NZ;z++) //fill voxels
    for (int y=0;y<NY;y++) {
      int zy = ((y+lo[1]) * inX) + ((z+lo[2]) * inXY);
      for (int x=0;x<NX;x++) {
        int j = lo[0] + x + zy;
        im[i] = img[j] - isolevel;
        i++;
      }
    }
  set_method(mcp, originalMC );
  run(mcp) ;
  clean_temps(mcp) ;
  if ((mcp->nverts < 3) || (mcp->ntrigs < 1)) {
    clean_all(mcp);
    free(mcp);
    return EXIT_FAILURE;
  }
  int npt = mcp->nverts;
  *vs = malloc(npt*sizeof(vec3d));
  for (int i = 0; i < npt; i++) {
    (*vs)[i].x = mcp->vertices[i].x + lo[0];
    (*vs)[i].y = mcp->vertices[i].y + lo[1];
    (*vs)[i].z = mcp->vertices[i].z + lo[2];
  }
  int ntri = mcp->ntrigs;
  *ts = malloc(ntri * sizeof(vec3i));
  for (int i=0;i<ntri;i++) {
    (*ts)[i].x = mcp->triangles[i].v3;
    (*ts)[i].y = mcp->triangles[i].v2;
    (*ts)[i].z = mcp->triangles[i].v1;
  }
  *nv = npt; //number of vertices
  *nt = ntri; //number of triangles
  clean_all(mcp);
  free(mcp);
  return EXIT_SUCCESS;
}
#endif //#ifdef NII2MESH

#ifdef MC_SELF_TEST
//These functions are only used by the self testing executable

static int littleEndianPlatform () {
  uint32_t value = 1;
  return (*((char *) &value) == 1);
}

void writePLY( MCB *mcb , const char *fn) {
  typedef struct  __attribute__((__packed__)) {
    uint8_t n;
    int32_t x,y,z;
  } vec1b3i;
  typedef struct {
    float x,y,z;
  } vec3s; //single precision (float32)
  int npt = mcb->nverts;
  int ntri = mcb->ntrigs;
  if ((npt < 3) || (ntri < 1)) {
    printf("Unable to create PLY file: No geometry\n");
    return;
  }
  FILE *fp = fopen(fn,"wb");
  if (fp == NULL)
    return;// EXIT_FAILURE;
  fputs("ply\n",fp);
  if (&littleEndianPlatform)
    fputs("format binary_little_endian 1.0\n",fp);
  else
    fputs("format binary_big_endian 1.0\n",fp);
  fputs("comment niimath\n",fp);
  char vpts[80];
  sprintf(vpts, "element vertex %d\n", npt);
  fwrite(vpts, strlen(vpts), 1, fp);
  fputs("property float x\n",fp);
  fputs("property float y\n",fp);
  fputs("property float z\n",fp);
  char vfc[80];
  sprintf(vfc, "element face %d\n", ntri);
  fwrite(vfc, strlen(vfc), 1, fp);
  fputs("property list uchar int vertex_indices\n",fp);
  fputs("end_header\n",fp);
  vec3s *pts32 = (vec3s *) malloc(npt * sizeof(vec3s));
  for (int i = 0; i < npt; i++) {
    pts32[i].x = mcb->vertices[i].x;
    pts32[i].y = mcb->vertices[i].y;
    pts32[i].z = mcb->vertices[i].z;
  }
  fwrite(pts32, npt * sizeof(vec3s), 1, fp);
  free(pts32);
  vec1b3i *tris4 = (vec1b3i *) malloc(ntri * sizeof(vec1b3i));
  for (int i = 0; i < ntri; i++) {
    tris4[i].n = 3;
    tris4[i].x = mcb->triangles[i].v3;
    tris4[i].y = mcb->triangles[i].v2;
    tris4[i].z = mcb->triangles[i].v1;
  }
  fwrite(tris4, ntri * sizeof(vec1b3i), 1, fp);
  free(tris4);
  fclose(fp);
}

typedef struct {
	int sizeof_hdr;
	char ignore[36];
	short dim[8];
	char ignore2[14];
	short datatype, bitpix, slice_start;
	float pixdim[8], vox_offset, scl_slope, scl_inter;
	char ignore3[224];
	char magic[4];
} TNIFTI;

float * readNIFTI(const char *fn, int *dim) {
  FILE* fp = fopen (fn, "r");
  TNIFTI hdr;
  size_t sz = fread(&hdr, sizeof(hdr), 1, fp);
  if ((hdr.sizeof_hdr != 348) || (hdr.vox_offset < 352)) {
    fclose(fp);
    printf("Error: only able to read uncompressed .nii NIfTI images of native endian.\n");
    return NULL;
  }
  dim[0] = hdr.dim[1];
  dim[1] = hdr.dim[2];
  dim[2] = hdr.dim[3];
  int datatype = hdr.datatype;
  int nvox = dim[0] * dim[1] * dim[2];
  fseek(fp, round(hdr.vox_offset), SEEK_SET);
  float * img32 = (float *) malloc(nvox*sizeof(float));
  if (datatype == 16) {
    size_t sz = fread(img32, nvox*4, 1, fp);
  } else if (datatype == 2) {
    uint8_t * img8 = (uint8_t *) malloc(nvox);
    size_t sz = fread(img8, nvox, 1, fp);
    for (int i = 0; i < nvox; i++)
      img32[i] = img8[i];
    free(img8);
  } else {
    fclose(fp);
    printf("Error: only able to read uint8 and float32 NIfTI images, not %d\n", datatype);
    return NULL;
  }
  fclose(fp);
  return img32;
}

void writeNIFTI(float *img32, int *dim, const char *fn ) {
  TNIFTI hdr;
  memset(&hdr, 0, sizeof(hdr));
  hdr.sizeof_hdr = 348;
  hdr.dim[0] = 3;
  hdr.dim[1] = dim[0];
  hdr.dim[2] = dim[1];
  hdr.dim[3] = dim[2];
  hdr.datatype = 16; //DT_FLOAT32
  hdr.bitpix = 32;
  hdr.pixdim[1] = 1.0;
  hdr.pixdim[2] = 1.0;
  hdr.pixdim[3] = 1.0;
  hdr.vox_offset = 352;
  hdr.scl_slope = 1.0;
  hdr.scl_inter = 0.0;
  hdr.magic[0] = 'n';
  hdr.magic[1] = '+';
  hdr.magic[2] = '1';
  FILE *fp = fopen( fn, "wb" ) ;
  fwrite(&hdr, sizeof(TNIFTI), 1, fp);
  int32_t dummy = 0;
  fwrite(&dummy, sizeof(dummy), 1, fp);
  int nvox = dim[0] * dim[1] * dim[2];
  fwrite(img32, nvox * sizeof(float), 1, fp);
  fclose(fp);
}

//_____________________________________________________________________________
// Compute data
float * compute_data(int obj_type, int *dim) {
  dim[0] = 60;
  dim[1] = 60;
  dim[2] = 60;
  float x,y,z      ;
  float sx,sy,sz   ;
  float tx,ty,tz , val = 0 ;
  int i, j, k;
  float a, b, c, d;
  float r = 1.85f ;
  float R = 4 ;
  if (obj_type < 0 || obj_type > 9) {
   fprintf(stderr,"Bad obj_type. Value must be between 0 and 9\n");
   return NULL;
  }
  sx     = (float) dim[0] / 16 ;
  sy     = (float) dim[1] / 16 ;
  sz     = (float) dim[2] / 16 ;
  tx     = (float) dim[0] / (2*sx) ;
  ty     = (float) dim[1] / (2*sy) + 1.5f ;
  tz     = (float) dim[2] / (2*sz) ;
  int vox = 0;
  float * img32 = (float *) malloc(dim[0]*dim[1]*dim[2]*sizeof(float));
  for(  k = 0 ; k < dim[2] ; k++ )
  {
    z = ( (float) k ) / sz  - tz ;
    for(  j = 0 ; j < dim[1] ; j++ )
    {
      y = ( (float) j ) / sy  - ty ;
      for(  i = 0 ; i < dim[0] ; i++ )
      {
        x = ( (float) i ) / sx - tx ;
        switch(obj_type) {
         case 0: //cushin
            val = z*z*x*x - z*z*z*z - 2*z*x*x + 2*z*z*z + x*x - z*z - (x*x - z)*(x*x - z) - y*y*y*y - 2*x*x*y*y - y*y*z*z + 2*y*y*z + y*y;
            break;
         case 1: //sphere
            val = ( (x-2)*(x-2) + (y-2)*(y-2) + (z-2)*(z-2) - 1 ) * ( (x+2)*(x+2) + (y-2)*(y-2) + (z-2)*(z-2) - 1 ) * ( (x-2)*(x-2) + (y+2)*(y+2) + (z-2)*(z-2) - 1 );
            break;
         case 2: // plane
            val = x+y+z -3;
            break;
         case 3: // cassini
            val = (x*x + y*y + z*z + 0.45f*0.45f)*(x*x + y*y + z*z + 0.45f*0.45f) - 16*0.45f*0.45f*(x*x + z*z) - 0.5f*0.5f;
            break;
         case 4: // blooby
            val = x*x*x*x - 5*x*x+ y*y*y*y - 5*y*y + z*z*z*z - 5*z*z + 11.8;
            break;
         case 5: // chair
            val = (x*x+y*y+z*z-0.95f*25)*(x*x+y*y+z*z-0.95f*25)-0.8f*((z-5)*(z-5)-2*x*x)*((z+5)*(z+5)-2*y*y);
            break;
         case 6: // cyclide
            b = 2; d = 6; a = 2; c = 3;
            val = ( x*x + y*y + z*z + b*b - d*d ) * ( x*x + y*y + z*z + b*b - d*d ) - 4 * ( ( a*x - c*d ) * ( a*x - c*d ) + b*b * y*y );
            break;
         case 7: // 2 torus
            val = ( ( x*x + y*y + z*z + R*R - r*r ) * ( x*x + y*y + z*z + R*R - r*r ) - 4 * R*R * ( x*x + y*y ) ) *
                  ( ( x*x + (y+R)*(y+R) + z*z + R*R - r*r ) * ( x*x + (y+R)*(y+R) + z*z + R*R - r*r ) - 4 * R*R * ( (y+R)*(y+R) + z*z ) ) ;
            break;
         case 8: // mc case
            val =    - 26.5298*(1-x)*(1-y)*(1-z) + 81.9199*x*(1-y)*(1-z) - 100.68*x*y*(1-z) + 3.5498*(1-x)*y*(1-z)
                     + 24.1201*(1-x)*(1-y)*  z   - 74.4702*x*(1-y)*  z   + 91.5298*x*y*  z  - 3.22998*(1-x)*y*  z ;
            break;
         case 9: // Drip
            val = x*x + y*y - 0.5*( 0.995*z*z + 0.005 - z*z*z ) +0.0025;  // -0.0754+0.01, -0.0025 + 0.01, grid 40^3, [-1.5,1.5];
            break;
        }
        img32[vox] = val;
        vox++;
      }
    }
  }
  #define MAX_BUF 1024
  char path[MAX_BUF], fn[MAX_BUF];
  if (getcwd (path, MAX_BUF) != path) return img32;
  int len = snprintf (fn, MAX_BUF-1, "%s/%d.nii", path, obj_type);
  if (len < 0) exit(EXIT_FAILURE);
  printf(  "Assuming write permissions to Save object %d as NIFTI %s\n", obj_type, fn);
  writeNIFTI(img32, dim, fn);
  return img32;
}

double clockMsec() { //return milliseconds since midnight
  struct timespec _t;
  clock_gettime(CLOCK_MONOTONIC, &_t);
  return _t.tv_sec*1000.0 + (_t.tv_nsec/1.0e6);
}

long timediff(double startTimeMsec, double endTimeMsec) {
  return round(endTimeMsec - startTimeMsec);
}

void self_test(const char *fn, float isolevel, int originalMC) {
  int dim[3];
  float * img32 = NULL;
  int obj_type = round(isolevel);
  if (strlen(fn) > 0) {
    img32 = readNIFTI(fn, dim);
    if (img32 == NULL) return;
    int nvox = dim[0]*dim[1]*dim[2];
    for (int i = 0; i < nvox; i++)
      img32[i] = img32[i] - isolevel;
  } else
    img32 = compute_data(obj_type, dim);
  if (img32 == NULL) return;
  MCB *mcp;
  mcp = MarchingCubes(-1, -1, -1);
  set_resolution( mcp, dim[0], dim[1], dim[2]) ;
  init_all(mcp) ;
  printf("input voxels: %d %d %d\n", dim[0], dim[1], dim[2]);
  set_method(mcp, originalMC );
  memcpy(mcp->data, img32, dim[0]*dim[1]*dim[2]*sizeof(float) ) ;
  free(img32);
  double startTime = clockMsec();
  run(mcp);
  printf("output mesh vert: %d tris: %d ms: %ld\n", mcp->nverts, mcp->ntrigs,  timediff(startTime, clockMsec()));
  clean_temps(mcp) ;
  #define MAX_BUF 1024
  char path[MAX_BUF], meshfn[MAX_BUF];
  if (getcwd (path, MAX_BUF) != path) exit(EXIT_FAILURE);
  int len = snprintf (meshfn, MAX_BUF-1, "%s/%d.ply", path, obj_type);
  if (len < 0) exit(EXIT_FAILURE);
  writePLY(mcp, meshfn) ;
  clean_all(mcp) ;
  free(mcp);
}

int main(int argc, char **argv) {
  if (argc > 2) {
    float isolevel = atof(argv[2]);
    int originalMC = 0;
    if (argc > 3)
      originalMC = atoi(argv[3]);
    self_test(argv[1], isolevel, originalMC);
    exit(EXIT_SUCCESS);
  }
  printf("Create PLY-format mesh from NIfTI format volume\n");
  printf("Usage:\n");
  printf("  %s niiname isolevel originalMC\n", argv[0]);
  printf("    niiname: uncompressed NIfTI .nii format volume\n");
  printf("    isolevel: boundary level for surface\n");
  printf("    originalMC: use Lewiner (0, default) or original (1) method \n");
  printf("Example:\n");
  printf("  %s bet.nii 64 1\n", argv[0]);
  printf("Running with no arguments creates test data\n");
  for (int i = 0; i < 10; i++)
    self_test("", i, 0);
  exit(EXIT_SUCCESS);
}
#endif //MC-SELF-TEST
