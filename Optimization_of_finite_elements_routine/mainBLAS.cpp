#include "Connectivity.hpp"
#include "Gmsh.hpp"
#include "Nodes.hpp"
#include "Numerics.hpp"
#include "Tools.hpp"
#include <mkl.h> // for the clock
#include <time.h>
#include <omp.h>

// dans codes-project :
// make blas ; cd benchmark/cube/ ; gmsh mesh.geo -3 ; ../../dg ; cd ../..

/*void afficherVecteur(const std::vector<double>& vecteur) {
    std::cout << "[ ";
    for (const auto& element : vecteur) {
        std::cout << element << " ";
    }
    std::cout << "]" << std::endl;
}*/

int main(int argc, char **argv) {

  // ======================================================================
  // 1A] Set simulation parameters
  // ======================================================================

  int iOutputGmsh;               // Number of time steps between gmsh output
  int N;                         // Degree of polynomial bases
  double FinalTime;              // Final time of the simulation
  vector<int> _CoefDataBaseNum;  // Database for parameters
  vector<double> _CoefDataBase1; // Database for parameters
  vector<double> _CoefDataBase2; // Database for parameters

  // Read setup file
  ifstream setupfile("setup");
  if (!setupfile.is_open()) {
    cerr << "ERROR: Setup file not opened.\n";
    exit(1);
  }
  string line;
  while (setupfile) {
    setupfile >> line;
    if (line == "$OutputStep")
      setupfile >> iOutputGmsh;
    if (line == "$PolynomialDegree")
      setupfile >> N;
    if (line == "$Duration")
      setupfile >> FinalTime;
    if (line == "$Param") {
      int PARAM;
      setupfile >> PARAM;
      _CoefDataBaseNum.resize(PARAM);
      _CoefDataBase1.resize(PARAM);
      _CoefDataBase2.resize(PARAM);
      for (int i = 0; i < PARAM; i++) {
        setupfile >> _CoefDataBaseNum[i];
        setupfile >> _CoefDataBase1[i];
        setupfile >> _CoefDataBase2[i];
      }
    }
  }

  // ======================================================================
  // 1B] Load mesh
  // ======================================================================

  int K;              // Number of elements in the mesh
  vector<double> _VX; // Coordinate 'x' of the vertices of the mesh [#vert]
  vector<double> _VY; // Coordinate 'y' of the vertices of the mesh [#vert]
  vector<double> _VZ; // Coordinate 'z' of the vertices of the mesh [#vert]
  vector<int> _EMsh;  // List of orginal gmsh numbers of elements [K]
  vector<int> _ETag;  // List of element tag (gmsh: physical tag) [K]
  vector<int> _EToV;  // Element-to-vertex connectivity array [K,NVertTet]

  // Load mesh number i for polynomial solution x^i
  char meshN[20]; // Définir une taille suffisante pour contenir le nom de fichier
  sprintf(meshN, "mesh.msh"); 
  // pour un maillage adapté au degré polynomial : sprintf(meshN, "mesh%1d.msh", N);
  // sinon sprintf(meshN, "mesh.msh");
  loadMeshGmsh(meshN, K, _VX, _VY, _VZ, _EMsh, _ETag, _EToV);

  // ======================================================================
  // 1C] Get finite element nodes
  // ======================================================================

  int Np = (N + 1) * (N + 2) * (N + 3) / 6; // Number of nodes per element
  int Nfp = (N + 1) * (N + 2) / 2;          // Number of nodes per face
  vector<double> _r, _s, _t;                // Coordinates of reference nodes [Np]
  vector<double> _x, _y, _z;                // Coordinates of physical nodes [K,Np]

  // Get local coordinates of nodes on the reference element
  getReferenceNodes(N, _r, _s, _t);

  // Get global coordinates of nodes on the physical elements
  getPhysicalNodes(_EToV, _VX, _VY, _VZ, _r, _s, _t, _x, _y, _z);

  // ======================================================================
  // 1D] Build connectivity matrices
  // ======================================================================

  vector<int> _EToE;  // Element-to-Element [K,NFacesTet]
  vector<int> _EToF;  // Element-to-Face [K,NFacesTet]
  vector<int> _NfToN; // LocalFaceNode-to-LocalNode [NFacesTet,Nfp]
  vector<int> _mapP;  // GlobalFaceNode-to-NeighborGlobalNode [K,Nfp,NFacesTet]

  // Build connectivity matrices
  buildConnectivity(K, Nfp, Np, _r, _s, _t, _x, _y, _z, _EToV, _EToE, _EToF, _NfToN, _mapP);

  // ======================================================================
  // 1E] Build geometric factors & elemental matrices
  // ======================================================================

  vector<double> _rstxyz; // ... for volume terms  [K,9]
  vector<double> _Fscale; // ... for surface terms [K,NFacesTet]
  vector<double> _nx;     // x-component of normal to the face [K,NFacesTet]
  vector<double> _ny;     // y-component of normal to the face [K,NFacesTet]
  vector<double> _nz;     // z-component of normal to the face [K,NFacesTet]
  vector<double> _Dr;     // r-derivative matrix [Np,Np]
  vector<double> _Ds;     // s-derivative matrix [Np,Np]
  vector<double> _Dt;     // t-derivative matrix [Np,Np]
  vector<double> _LIFT;   // Lift matrix [NFacesTet,Nfp,Np]

  // Build geometric factors
  buildGeomFactors(_EToV, _VX, _VY, _VZ, _rstxyz, _Fscale, _nx, _ny, _nz);

  // Build elemental matrices
  buildElemMatrices(N, _r, _s, _t, _NfToN, _Dr, _Ds, _Dt, _LIFT);

  // ======================================================================
  // 1F] Build physical coefficient maps
  // ======================================================================

  vector<double> _c(K);
  vector<double> _rho(K);

  for (int k = 0; k < K; k++) {
    int i = 0;
    while (i < _CoefDataBaseNum.size()) {
      if (_CoefDataBaseNum[i] == _ETag[k]) {
        _c[k] = _CoefDataBase1[i];
        _rho[k] = _CoefDataBase2[i];
        break;
      }
      i++;
    }
    if (i == _CoefDataBaseNum.size()) {
      cerr << "ERROR: ETag " << _ETag[k] << " not found in parameter database\n";
      exit(1);
    }
  }
  exportParamGmsh("info_c", _EMsh, _c);
  exportParamGmsh("info_rho", _EMsh, _rho);

  // ======================================================================
  // 1G] Build time stepping scheme
  // ======================================================================

  double dt = 1e9;
  for (int k = 0; k < K; k++) {
    for (int f = 0; f < NFacesTet; ++f) {
      double val = 1. / (_c[k] * (N + 1) * (N + 1) * _Fscale[k * NFacesTet + f]);
      if (val < dt)
        dt = val;
    }
  }
  double CFL = 0.75;                 // CFL
  dt *= CFL;                         // Time step
  int Nsteps = ceil(FinalTime / dt); // Number of global time steps

  //afficherVecteur(_Fscale);
  //afficherVecteur(_c);

  vector<double> rk4a, rk4b, rk4c;
  getRungeKuttaCoefficients(rk4a, rk4b, rk4c);

  cout << "INFO: DT = " << dt << " and Nsteps = " << Nsteps << "\n";

  // ======================================================================
  // 1H] Memory storage for fields, RHS and residual at nodes
  // ======================================================================

  int Nfields = 4; // Number of unknown fields

  // Memory storage at each node of each element
  vector<double> _valQ(K * Np * Nfields, 0.); // Values of fields
  vector<double> _rhsQ(K * Np * Nfields, 0.); // RHS
  vector<double> _resQ(K * Np * Nfields, 0.); //

  // Initialization
  for (int k = 0; k < K; ++k) {
    for (int n = 0; n < Np; ++n) {
      double x = _x[k * Np + n];
      double y = _y[k * Np + n];
      double z = _z[k * Np + n];
      _valQ[k * Np * Nfields + n * Nfields + 0] = exp(-(x * x + y * y + z * z) / 0.1);
    }
  }

  cout << "INFO: K * Np * Nfields = " << K * Np * Nfields << "\n";

  // Export initial solution
  if (iOutputGmsh > 0)
    exportSolGmsh(N, _r, _s, _t, _EMsh, 0, 0., _valQ);

  // =========================================================================================================
  // 2] RUN
  // =========================================================================================================

  // estimation de performance : départ boucle temporelle
  double timeBegin = dsecnd();

  // Global time iteration
  for (int nGlo = 0; nGlo < Nsteps; ++nGlo) {
    double runTime = nGlo * dt; // Time at the beginning of the step

    // Local time iteration
    for (int nLoc = 0; nLoc < 5; ++nLoc) {
      double a = rk4a[nLoc];
      double b = rk4b[nLoc];

      // ======================== (1) UPDATE RHS

      for (int k = 0; k < K; ++k) {

        double c = _c[k];
        double rho = _rho[k];

        // ======================== (1.1) UPDATE PENALTY VECTOR

        double* s_p_flux = (double *)(malloc(Nfp * NFacesTet * sizeof(double)));
        double* s_u_flux = (double *)(malloc(Nfp * NFacesTet * sizeof(double)));
        double* s_v_flux = (double *)(malloc(Nfp * NFacesTet * sizeof(double)));
        double* s_w_flux = (double *)(malloc(Nfp * NFacesTet * sizeof(double)));

        for (int f = 0; f < NFacesTet; f++) {

          // Fetch normal
          double nx = _nx[k * NFacesTet + f];
          double ny = _ny[k * NFacesTet + f];
          double nz = _nz[k * NFacesTet + f];

          // Fetch medium parameters
          int k2 = _EToE[k * NFacesTet + f];
          double cP = _c[k2];
          double rhoP = _rho[k2];

          // Compute penalty terms
          for (int nf = f * Nfp; nf < (f + 1) * Nfp; nf++) {

            // Index of node in current element
            int n1 = _NfToN[nf];

            // Index of node in neighbor element
            int n2 = _mapP[nf + k * NFacesTet * Nfp];

            // Load values 'minus' corresponding to current element
            double pM = _valQ[k * Np * Nfields + n1 * Nfields + 0];
            double uM = _valQ[k * Np * Nfields + n1 * Nfields + 1];
            double vM = _valQ[k * Np * Nfields + n1 * Nfields + 2];
            double wM = _valQ[k * Np * Nfields + n1 * Nfields + 3];
            double nMdotuM = (nx * uM + ny * vM + nz * wM);

            if (n2 >= 0) { // ... if there is a neighbor element ...

              // Load values 'plus' corresponding to neighbor element
              double pP = _valQ[k2 * Np * Nfields + n2 * Nfields + 0];
              double uP = _valQ[k2 * Np * Nfields + n2 * Nfields + 1];
              double vP = _valQ[k2 * Np * Nfields + n2 * Nfields + 2];
              double wP = _valQ[k2 * Np * Nfields + n2 * Nfields + 3];
              double nMdotuP = (nx * uP + ny * vP + nz * wP);

              // Penalty terms for interface between two elements
              double temp = c / (rhoP * cP + rho * c) * ((pP - pM) - rhoP * cP * (nMdotuP - nMdotuM));
              s_p_flux[nf] = c / (1. / (rhoP * cP) + 1. / (rho * c)) * ((nMdotuP - nMdotuM) - 1. / (rhoP * cP) * (pP - pM));
              s_u_flux[nf] = nx * temp;
              s_v_flux[nf] = ny * temp;
              s_w_flux[nf] = nz * temp;

            } else {

              // Homogeneous Dirichlet on 'p'
              double tmp = -2. / (rhoP * cP) * pM;
              // Homogeneous Dirichlet on 'u'
              // double tmp = 2*nMdotuM;
              // ABC
              // double tmp = nMdotuM - 1./(rhoP*cP) * pM;

              // Penalty terms for boundary of the domain
              s_p_flux[nf] = -c / (1. / (rhoP * cP) + 1. / (rho * c)) * tmp;
              s_u_flux[nf] = nx * c / (rhoP * cP + rho * c) * tmp;
              s_v_flux[nf] = ny * c / (rhoP * cP + rho * c) * tmp;
              s_w_flux[nf] = nz * c / (rhoP * cP + rho * c) * tmp;
            }
          }
        }

        // ======================== (1.2) COMPUTING VOLUME TERMS

        // Load geometric factors
        double rx = _rstxyz[k * 9 + 0];
        double ry = _rstxyz[k * 9 + 1];
        double rz = _rstxyz[k * 9 + 2];
        double sx = _rstxyz[k * 9 + 3];
        double sy = _rstxyz[k * 9 + 4];
        double sz = _rstxyz[k * 9 + 5];
        double tx = _rstxyz[k * 9 + 6];
        double ty = _rstxyz[k * 9 + 7];
        double tz = _rstxyz[k * 9 + 8];

        // Load fields
        double* s_p = (double *)(malloc(Np * sizeof(double)));;
        double* s_u = (double *)(malloc(Np * sizeof(double)));;
        double* s_v = (double *)(malloc(Np * sizeof(double)));;
        double* s_w = (double *)(malloc(Np * sizeof(double)));;
        for (int n = 0; n < Np; ++n) {
          s_p[n] = _valQ[k * Np * Nfields + n * Nfields + 0];
          s_u[n] = _valQ[k * Np * Nfields + n * Nfields + 1];
          s_v[n] = _valQ[k * Np * Nfields + n * Nfields + 2];
          s_w[n] = _valQ[k * Np * Nfields + n * Nfields + 3];
        }

        // Compute mat-vec product for surface term
        double* dpdr = (double *)(malloc(Np * sizeof(double)));;
        double* dpds = (double *)(malloc(Np * sizeof(double)));;
        double* dpdt = (double *)(malloc(Np * sizeof(double)));;
        double* dudr = (double *)(malloc(Np * sizeof(double)));;
        double* duds = (double *)(malloc(Np * sizeof(double)));;
        double* dudt = (double *)(malloc(Np * sizeof(double)));;
        double* dvdr = (double *)(malloc(Np * sizeof(double)));;
        double* dvds = (double *)(malloc(Np * sizeof(double)));;
        double* dvdt = (double *)(malloc(Np * sizeof(double)));;
        double* dwdr = (double *)(malloc(Np * sizeof(double)));;
        double* dwds = (double *)(malloc(Np * sizeof(double)));;
        double* dwdt = (double *)(malloc(Np * sizeof(double)));;
        cblas_dgemv(CblasRowMajor, CblasTrans, Np, Np, 1.0, _Dr.data(), Np, s_p, 1, 1.0, dpdr, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, Np, Np, 1.0, _Ds.data(), Np, s_p, 1, 1.0, dpds, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, Np, Np, 1.0, _Dt.data(), Np, s_p, 1, 1.0, dpdt, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, Np, Np, 1.0, _Dr.data(), Np, s_u, 1, 1.0, dudr, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, Np, Np, 1.0, _Ds.data(), Np, s_u, 1, 1.0, duds, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, Np, Np, 1.0, _Dt.data(), Np, s_u, 1, 1.0, dudt, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, Np, Np, 1.0, _Dr.data(), Np, s_v, 1, 1.0, dvdr, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, Np, Np, 1.0, _Ds.data(), Np, s_v, 1, 1.0, dvds, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, Np, Np, 1.0, _Dt.data(), Np, s_v, 1, 1.0, dvdt, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, Np, Np, 1.0, _Dr.data(), Np, s_w, 1, 1.0, dwdr, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, Np, Np, 1.0, _Ds.data(), Np, s_w, 1, 1.0, dwds, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, Np, Np, 1.0, _Dt.data(), Np, s_w, 1, 1.0, dwdt, 1);
        
        double *dpdx = (double *)malloc(Np * sizeof(double));
        double *dpdy = (double *)malloc(Np * sizeof(double));
        double *dpdz = (double *)malloc(Np * sizeof(double));
        double *dudx = (double *)malloc(Np * sizeof(double));
        double *dvdy = (double *)malloc(Np * sizeof(double));
        double *dwdz = (double *)malloc(Np * sizeof(double));
        double *divU = (double *)malloc(Np * sizeof(double));

        cblas_daxpy(Np, rx, dpdr, 1, dpdx, 1);
        cblas_daxpy(Np, sx, dpds, 1, dpdx, 1);
        cblas_daxpy(Np, tx, dpdt, 1, dpdx, 1);

        cblas_daxpy(Np, ry, dpdr, 1, dpdy, 1);
        cblas_daxpy(Np, sy, dpds, 1, dpdy, 1);
        cblas_daxpy(Np, ty, dpdt, 1, dpdy, 1);

        cblas_daxpy(Np, rz, dpdr, 1, dpdz, 1);
        cblas_daxpy(Np, sz, dpds, 1, dpdz, 1);
        cblas_daxpy(Np, tz, dpdt, 1, dpdz, 1);

        cblas_daxpy(Np, rx, dudr, 1, dudx, 1);
        cblas_daxpy(Np, sx, duds, 1, dudx, 1);
        cblas_daxpy(Np, tx, dudt, 1, dudx, 1);

        cblas_daxpy(Np, ry, dvdr, 1, dvdy, 1);
        cblas_daxpy(Np, sy, dvds, 1, dvdy, 1);
        cblas_daxpy(Np, ty, dvdt, 1, dvdy, 1);

        cblas_daxpy(Np, rz, dwdr, 1, dwdz, 1);
        cblas_daxpy(Np, sz, dwds, 1, dwdz, 1);
        cblas_daxpy(Np, tz, dwdt, 1, dwdz, 1);

        cblas_dcopy(Np, dudx, 1, divU, 1);
        cblas_daxpy(Np, 1.0, dvdy, 1, divU, 1);
        cblas_daxpy(Np, 1.0, dwdz, 1, divU, 1);

        for (int n = 0; n < Np; ++n) {
          // Compute RHS (only part corresponding to volume terms)
          _rhsQ[k * Np * Nfields + n * Nfields + 0] = -c * c * rho * divU[n];
          _rhsQ[k * Np * Nfields + n * Nfields + 1] = -1. / rho * dpdx[n];
          _rhsQ[k * Np * Nfields + n * Nfields + 2] = -1. / rho * dpdy[n];
          _rhsQ[k * Np * Nfields + n * Nfields + 3] = -1. / rho * dpdz[n];
        }

        free(dpdr);
        free(dpds);
        free(dpdt);
        free(dudr);
        free(duds);
        free(dudt);
        free(dvdr);
        free(dvds);
        free(dvdt);
        free(dwdr);
        free(dwds);
        free(dwdt);
    
        free(s_p);
        free(s_u);
        free(s_v);
        free(s_w);

        // ======================== (1.3) COMPUTING SURFACE TERMS

          for (int f = 0; f < NFacesTet; f++) { // boucle sur les lignes de _rhsQ
            double* pU_lift = (double *)(malloc(4*Nfp*NFacesTet*sizeof(double)));
            double Fscale = _Fscale[k * NFacesTet + f];
            double* _rhsQ_f = &_rhsQ[k*Np*Nfields];
            double* _LIFT_f = &_LIFT[f*Nfp];

            for (int m=0 ; m<Nfp ; ++m){
              pU_lift[4*m+0] = s_p_flux[m];
              pU_lift[4*m+1] = s_u_flux[m];
              pU_lift[4*m+2] = s_v_flux[m];
              pU_lift[4*m+3] = s_w_flux[m];
            }
            // Compute mat-vec product for surface term
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,Np,4,Nfp,-Fscale,&_LIFT[f*Nfp],NFacesTet*Nfp,&pU_lift[Nfp*f],4,1.0,&_rhsQ[k*Np*Nfields],4);
            free(pU_lift);
          }
      }
      // ======================== (2) UPDATE RESIDUAL + FIELDS

      for (int k = 0; k < K; ++k) {
        for (int n = 0; n < Np; ++n) {
          for (int iField = 0; iField < Nfields; ++iField) {
            int id = k * Np * Nfields + n * Nfields + iField;
            _resQ[id] = a * _resQ[id] + dt * _rhsQ[id];
            _valQ[id] = _valQ[id] + b * _resQ[id];
          }
        }
      }
    }

    // Export solution
    if ((iOutputGmsh > 0) && ((nGlo + 1) % iOutputGmsh == 0))
      exportSolGmsh(N, _r, _s, _t, _EMsh, nGlo + 1, runTime + dt, _valQ);
  }

  // timing : end and calculus of performance
  double timeEnd = dsecnd();
  double duration = timeEnd - timeBegin;

  double Nloc_run = 5;
  double NFacesTet_run = 4; // nombre de faces des éléments de référence
  double nb_op = Nsteps*Nloc_run*K*(22 + NFacesTet_run*(8 + Nfp*80) + Np*(75 + Np*30 + NFacesTet_run*(8 + Nfp*8)));
  double debit_arith = nb_op / duration;
  cout << "debit : " << debit_arith*1e-9 << " GFLOPS ; temps de calcul : " << duration << " sec "<< endl;

  // ======================================================================
  // 3] Post-processing
  // ======================================================================

  // Export final solution
  if (iOutputGmsh > 0)
    exportSolGmsh(N, _r, _s, _t, _EMsh, Nsteps, Nsteps * dt, _valQ);

  return 0;
}
