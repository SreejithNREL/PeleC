#include "prob.H"

void
init_bc()
{
  amrex::Real vt, ek, a, yl, yr, sumY, T, rho, e;
  amrex::Real molefrac[NUM_SPECIES], massfrac[NUM_SPECIES];

  a = 2.0;
  for (int n = 0; n < NUM_SPECIES; n++)
    molefrac[n] = 0.0;

 molefrac[O2_ID]  = 1.0 / (1.0 + PeleC::h_prob_parm_device->phi_in / a + 0.79 / 0.21);
 molefrac[CH4_ID] = PeleC::h_prob_parm_device->phi_in * molefrac[O2_ID] / a;
 molefrac[N2_ID]  = 1.0 - molefrac[CH4_ID] - molefrac[O2_ID];

  // molefrac[N2_ID] = 0.79;
  // molefrac[O2_ID] = 0.21;
  T = PeleC::h_prob_parm_device->T_in;

  double p = PeleC::h_prob_parm_device->pamb;
  auto eos = pele::physics::PhysicsType::eos();

  eos.X2Y(molefrac, massfrac);
  eos.PYT2RE(p, massfrac, T, rho, e);

  vt = PeleC::h_prob_parm_device->vn_in;
  ek = 0.5 * (vt * vt);

  (PeleC::h_prob_parm_device->fuel_state)[URHO]  = rho;
  (PeleC::h_prob_parm_device->fuel_state)[UMX]   = rho * vt;
  (PeleC::h_prob_parm_device->fuel_state)[UMY]   = 0.0;
  (PeleC::h_prob_parm_device->fuel_state)[UMZ]   = 0.0;
  (PeleC::h_prob_parm_device->fuel_state)[UEINT] = rho * e;
  (PeleC::h_prob_parm_device->fuel_state)[UEDEN] = rho * (e + ek);
  (PeleC::h_prob_parm_device->fuel_state)[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++){
    (PeleC::h_prob_parm_device->fuel_state)[UFS + n] = rho * massfrac[n];
  }

  p = 990.0e+6;
  T = 25500.;
  eos.X2Y(molefrac, massfrac);
  eos.PYT2RE(p, massfrac, T, rho, e);
  (PeleC::h_prob_parm_device->kernel_state)[URHO]  = rho;
  (PeleC::h_prob_parm_device->kernel_state)[UMX]   = rho * vt;
  (PeleC::h_prob_parm_device->kernel_state)[UMY]   = 0.0;
  (PeleC::h_prob_parm_device->kernel_state)[UMZ]   = 0.0;
  (PeleC::h_prob_parm_device->kernel_state)[UEINT] = rho * e;
  (PeleC::h_prob_parm_device->kernel_state)[UEDEN] = rho * (e + ek);
  (PeleC::h_prob_parm_device->kernel_state)[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++)
    (PeleC::h_prob_parm_device->kernel_state)[UFS + n] = rho * massfrac[n];
}

void
pc_prob_close()
{

}

extern "C" 
{
  void amrex_probinit(
              const int* init,
              const int* name,
              const int* namelen,
              const amrex_real* problo,
              const amrex_real* probhi)
    {

      amrex::ParmParse pp("prob");
      pp.query("pamb"       , PeleC::h_prob_parm_device->pamb);
      pp.query("phi_in"     , PeleC::h_prob_parm_device->phi_in);
      pp.query("T_in"       , PeleC::h_prob_parm_device->T_in);
      pp.query("vn_in"      , PeleC::h_prob_parm_device->vn_in);
      pp.query("turbulence" , PeleC::h_prob_parm_device->turbulence);
      pp.query("init_kernel", PeleC::h_prob_parm_device->init_kernel);
      pp.query("pertmag"    , PeleC::h_prob_parm_device->pertmag);
      pp.query("iname"      , PeleC::prob_parm_host->iname);
      pp.query("turb_ic"    , PeleC::prob_parm_host->turb_ic);

      PeleC::h_prob_parm_device->L[0] = probhi[0] - problo[0];
      PeleC::h_prob_parm_device->L[1] = probhi[1] - problo[1];
      PeleC::h_prob_parm_device->L[2] = probhi[2] - problo[2];

      // PeleC::h_prob_parm_device->fuel_state.resize(NVAR);
      // PeleC::h_prob_parm_device->kernel_state.resize(NVAR);

      if (PeleC::h_prob_parm_device->restart) {
        amrex::Print() << "Skipping input file reading and assuming restart."
                       << std::endl;
      } 
      else if(PeleC::h_prob_parm_device->init_kernel == false){
        //initialize premixed state
        init_bc();
      }
      else {
        //get array sizes
        int z_coord = 0;
        get_inputs(PeleC::prob_parm_host->iname,z_coord);
        
        PeleC::h_prob_parm_device->nx = PeleC::prob_parm_host->nx_tmp;
        PeleC::h_prob_parm_device->ny = PeleC::prob_parm_host->ny_tmp;
        PeleC::h_prob_parm_device->nz = PeleC::prob_parm_host->nz_tmp;
        PeleC::h_prob_parm_device->nscal = PeleC::prob_parm_host->nscal_tmp;

        amrex::Vector<double> data_tmp(PeleC::h_prob_parm_device->nx * PeleC::h_prob_parm_device->ny * PeleC::h_prob_parm_device->nscal);
        amrex::Vector<double> gridx_input(PeleC::h_prob_parm_device->nx);
        amrex::Vector<double> gridy_input(PeleC::h_prob_parm_device->ny);
        amrex::Vector<double> gridz_input(PeleC::h_prob_parm_device->nz);

        // //get data from file
        read_input_file(PeleC::prob_parm_host->iname,
                        z_coord,
                        PeleC::h_prob_parm_device->nx,
                        PeleC::h_prob_parm_device->ny,
                        PeleC::h_prob_parm_device->nz,
                        PeleC::h_prob_parm_device->nscal,
                        gridx_input, gridy_input, gridz_input, data_tmp);

        //copy data to the working arrays
        PeleC::prob_parm_host->v_xinput.resize(PeleC::h_prob_parm_device->nx);
        for (int i = 0; i < PeleC::h_prob_parm_device->nx; i++) {
          PeleC::prob_parm_host->v_xinput[i] = gridx_input[i];
        }

        PeleC::prob_parm_host->v_yinput.resize(PeleC::h_prob_parm_device->ny);
        for (int i = 0; i < PeleC::h_prob_parm_device->ny; i++) {
          PeleC::prob_parm_host->v_yinput[i] = gridy_input[i];
        }

        PeleC::prob_parm_host->data.resize(PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny*PeleC::h_prob_parm_device->nscal);
        for (int i = 0; i < PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny*PeleC::h_prob_parm_device->nscal; i++) {
          PeleC::prob_parm_host->data[i] = data_tmp[i];
          // printf("data_tmp[i] %f \n", data_tmp[i]);
        }


        // Get the array differences.
        PeleC::prob_parm_host->v_xdiff.resize(PeleC::h_prob_parm_device->nx);
        std::adjacent_difference(
          PeleC::prob_parm_host->v_xinput.begin(), PeleC::prob_parm_host->v_xinput.end(),
          PeleC::prob_parm_host->v_xdiff.begin());
        PeleC::prob_parm_host->v_xdiff[0] = PeleC::prob_parm_host->v_xdiff[1];

        PeleC::prob_parm_host->v_ydiff.resize(PeleC::h_prob_parm_device->ny);
        std::adjacent_difference(
          PeleC::prob_parm_host->v_yinput.begin(), PeleC::prob_parm_host->v_yinput.end(),
          PeleC::prob_parm_host->v_ydiff.begin());
        PeleC::prob_parm_host->v_ydiff[0] = PeleC::prob_parm_host->v_ydiff[1];

        // Get pointer to the data
        PeleC::prob_parm_host->xarray.resize(PeleC::prob_parm_host->v_xinput.size());
        PeleC::prob_parm_host->xdiff.resize(PeleC::prob_parm_host->v_xdiff.size());
        PeleC::prob_parm_host->yarray.resize(PeleC::prob_parm_host->v_yinput.size());
        PeleC::prob_parm_host->ydiff.resize(PeleC::prob_parm_host->v_ydiff.size());
        PeleC::prob_parm_host->data_ic.resize(PeleC::prob_parm_host->data.size());

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->v_xinput.begin(),
         PeleC::prob_parm_host->v_xinput.end(),
         PeleC::prob_parm_host->xarray.begin());
            
        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->v_xdiff.begin(),
         PeleC::prob_parm_host->v_xdiff.end(),
         PeleC::prob_parm_host->xdiff.begin());
            
        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->v_yinput.begin(),
         PeleC::prob_parm_host->v_yinput.end(),
         PeleC::prob_parm_host->yarray.begin());
            
        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->v_ydiff.begin(),
         PeleC::prob_parm_host->v_ydiff.end(),
         PeleC::prob_parm_host->ydiff.begin());
            
        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->data.begin(),
         PeleC::prob_parm_host->data.end(),
         PeleC::prob_parm_host->data_ic.begin());

        PeleC::h_prob_parm_device->d_xarray  = PeleC::prob_parm_host->xarray.data();
        PeleC::h_prob_parm_device->d_xdiff   = PeleC::prob_parm_host->xdiff.data();
        PeleC::h_prob_parm_device->d_yarray  = PeleC::prob_parm_host->yarray.data();
        PeleC::h_prob_parm_device->d_ydiff   = PeleC::prob_parm_host->ydiff.data();
        PeleC::h_prob_parm_device->d_data_ic = PeleC::prob_parm_host->data_ic.data();


        // Dimensions of the input box.
        PeleC::h_prob_parm_device->Lxinput =
          PeleC::prob_parm_host->v_xinput[PeleC::h_prob_parm_device->nx - 1] + 0.5 * PeleC::prob_parm_host->v_xdiff[PeleC::h_prob_parm_device->nx - 1];
        PeleC::h_prob_parm_device->Lyinput =
          PeleC::prob_parm_host->v_yinput[PeleC::h_prob_parm_device->ny - 1] + 0.5 * PeleC::prob_parm_host->v_ydiff[PeleC::h_prob_parm_device->ny - 1];

        //initialize premixed state
        init_bc();
      }

      if(PeleC::h_prob_parm_device->turbulence){
        amrex::Print() << "Initializing turbulence from: " <<  PeleC::prob_parm_host->turb_ic.c_str() << std::endl;

        //get array sizes
        int z_coord = 1;
        get_inputs(PeleC::prob_parm_host->turb_ic,z_coord);
        
        //define local variable containing size info
        int nx = PeleC::prob_parm_host->nx_tmp;
        int ny = PeleC::prob_parm_host->ny_tmp;
        int nz = PeleC::prob_parm_host->nz_tmp;
        int nscal = PeleC::prob_parm_host->nscal_tmp; 

        //initialize variables for interpolation later on
        PeleC::h_prob_parm_device->ires_x = nx;
        PeleC::h_prob_parm_device->ires_y = ny;
        PeleC::h_prob_parm_device->ires_z = nz;

        amrex::Vector<double> data_tmp(nx * ny * nz * nscal);

        amrex::Vector<double> gridx_input(nx);
        amrex::Vector<double> gridy_input(ny);
        amrex::Vector<double> gridz_input(nz);

        // //get data from file
        read_input_file(PeleC::prob_parm_host->turb_ic,
                        z_coord,
                        nx,
                        ny,
                        nz,
                        nscal,
                        gridx_input, gridy_input, gridz_input, data_tmp);
        
        //copy data to the working arrays
        PeleC::prob_parm_host->T_xinput.resize(PeleC::h_prob_parm_device->ires_x);
        for (int i = 0; i < PeleC::h_prob_parm_device->ires_x; i++) {
          PeleC::prob_parm_host->T_xinput[i] = gridx_input[i];
        }
        PeleC::prob_parm_host->T_yinput.resize(PeleC::h_prob_parm_device->ires_y);
        for (int i = 0; i < PeleC::h_prob_parm_device->ires_y; i++) {
          PeleC::prob_parm_host->T_yinput[i] = gridy_input[i];
        }
        PeleC::prob_parm_host->T_zinput.resize(PeleC::h_prob_parm_device->ires_z);
        for (int i = 0; i < PeleC::h_prob_parm_device->ires_z; i++) {
          PeleC::prob_parm_host->T_zinput[i] = gridz_input[i];
        }
        PeleC::prob_parm_host->data_turb.resize(nx*ny*nz*nscal);
        for (int i = 0; i < nx*ny*nz*nscal; i++) {
          PeleC::prob_parm_host->data_turb[i] = data_tmp[i];
        }


        // Get the array differences.
        PeleC::prob_parm_host->T_xdiff.resize(PeleC::h_prob_parm_device->ires_x);
        std::adjacent_difference(
          PeleC::prob_parm_host->T_xinput.begin(), PeleC::prob_parm_host->T_xinput.end(),
          PeleC::prob_parm_host->T_xdiff.begin());
        PeleC::prob_parm_host->T_xdiff[0] = PeleC::prob_parm_host->T_xdiff[1];

        PeleC::prob_parm_host->T_ydiff.resize(PeleC::h_prob_parm_device->ires_y);
        std::adjacent_difference(
          PeleC::prob_parm_host->T_yinput.begin(), PeleC::prob_parm_host->T_yinput.end(),
          PeleC::prob_parm_host->T_ydiff.begin());
        PeleC::prob_parm_host->T_ydiff[0] = PeleC::prob_parm_host->T_ydiff[1];

        PeleC::prob_parm_host->T_zdiff.resize(PeleC::h_prob_parm_device->ires_z);
        std::adjacent_difference(
          PeleC::prob_parm_host->T_zinput.begin(), PeleC::prob_parm_host->T_zinput.end(),
          PeleC::prob_parm_host->T_zdiff.begin());
        PeleC::prob_parm_host->T_zdiff[0] = PeleC::prob_parm_host->T_zdiff[1];


        PeleC::prob_parm_host->data_ic.resize(PeleC::prob_parm_host->data.size());

        PeleC::prob_parm_host->xarray_turb.resize(PeleC::prob_parm_host->T_xinput.size());
        PeleC::prob_parm_host->xdiff_turb.resize(PeleC::prob_parm_host->T_xdiff.size());
        PeleC::prob_parm_host->yarray_turb.resize(PeleC::prob_parm_host->T_yinput.size());
        PeleC::prob_parm_host->ydiff_turb.resize(PeleC::prob_parm_host->T_ydiff.size());
        PeleC::prob_parm_host->zarray_turb.resize(PeleC::prob_parm_host->T_zinput.size());
        PeleC::prob_parm_host->zdiff_turb.resize(PeleC::prob_parm_host->T_zdiff.size());
        PeleC::prob_parm_host->data_ic_turb.resize(PeleC::prob_parm_host->data_turb.size());


        // Get pointer to the data
        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->T_xinput.begin(),
        PeleC::prob_parm_host->T_xinput.end(),  PeleC::prob_parm_host->xarray_turb.begin()
          );

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->T_xdiff.begin(),
        PeleC::prob_parm_host->T_xdiff.end(), PeleC::prob_parm_host->xdiff_turb.begin()
          );

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->T_yinput.begin(),
        PeleC::prob_parm_host->T_yinput.end(),  PeleC::prob_parm_host->yarray_turb.begin()
          );

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->T_ydiff.begin(),
        PeleC::prob_parm_host->T_ydiff.end(), PeleC::prob_parm_host->ydiff_turb.begin()
          );

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->T_zinput.begin(),
        PeleC::prob_parm_host->T_zinput.end(),  PeleC::prob_parm_host->zarray_turb.begin()
          );

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->T_zdiff.begin(),
        PeleC::prob_parm_host->T_zdiff.end(), PeleC::prob_parm_host->zdiff_turb.begin()
          );

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->data_turb.begin(),
        PeleC::prob_parm_host->data_turb.end(), PeleC::prob_parm_host->data_ic_turb.begin()
          );


        PeleC::h_prob_parm_device->d_xarray_turb  = PeleC::prob_parm_host->xarray_turb.data();
        PeleC::h_prob_parm_device->d_xdiff_turb   = PeleC::prob_parm_host->xdiff_turb.data();
        PeleC::h_prob_parm_device->d_yarray_turb  = PeleC::prob_parm_host->yarray_turb.data();
        PeleC::h_prob_parm_device->d_ydiff_turb   = PeleC::prob_parm_host->ydiff_turb.data();
        PeleC::h_prob_parm_device->d_zarray_turb  = PeleC::prob_parm_host->zarray_turb.data();
        PeleC::h_prob_parm_device->d_zdiff_turb   = PeleC::prob_parm_host->zdiff_turb.data();
        PeleC::h_prob_parm_device->d_data_ic_turb = PeleC::prob_parm_host->data_ic_turb.data();


        // Dimensions of the input box
        PeleC::h_prob_parm_device->Lxturb =
          PeleC::prob_parm_host->T_xinput[nx - 1] + 0.5 * PeleC::prob_parm_host->T_xinput[nx - 1];

        PeleC::h_prob_parm_device->Lyturb =
          PeleC::prob_parm_host->T_yinput[ny - 1] + 0.5 * PeleC::prob_parm_host->T_yinput[ny - 1];

        PeleC::h_prob_parm_device->Lzturb =
          PeleC::prob_parm_host->T_zinput[nz - 1] + 0.5 * PeleC::prob_parm_host->T_zinput[nz - 1];
      }
    }
}

void
read_input_file(
  const std::string iname,
  int z_coord,
  int nx,
  int ny,
  int nz,
  int nscal,
  amrex::Vector<amrex::Real>& gridx_input,
  amrex::Vector<amrex::Real>& gridy_input,
  amrex::Vector<amrex::Real>& gridz_input,
  amrex::Vector<amrex::Real>& data)
{
  std::ifstream infile(iname, std::ios::in);
  const std::string memfile = read_file(infile);
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }
  infile.close();
  int i,j,k,n;
  
  std::istringstream iss(memfile);

  // Read the file
  int nlines = 0;
  std::string firstline, line;
  if(z_coord){
    for (i=0;i<nscal+7;i++){
      std::getline(iss, firstline); // skip header
    }  
  }
  else{
    for (i=0;i<nscal+5;i++){
      std::getline(iss, firstline); // skip header
    }    
  }


  // Read grid in X
  int cnt = 0;
  for (i=0;i<nx;i++){
    std::getline(iss, line);
    std::istringstream sinput(line);
    sinput >> gridx_input[cnt];

    if(z_coord == 0){
      gridx_input[cnt] = gridx_input[cnt]*100.0; //convert to [cm]
    }
    cnt++;
  }

  // Read the data from the file
  cnt = 0;
  for (j=0;j<ny;j++){
    std::getline(iss, line);
    std::istringstream sinput(line);
    sinput >> gridy_input[cnt];
    if(z_coord == 0){
      gridy_input[cnt] = gridy_input[cnt]*100.0; //convert to [cm]
    }
    cnt++;
  }
  
  if(z_coord){
    // Read the data from the file
    cnt = 0;
    for (k=0;k<nz;k++){
      std::getline(iss, line);
      std::istringstream sinput(line);
      sinput >> gridz_input[cnt];
      // gridz_input[cnt] = gridz_input[cnt]; //convert to [cm]
      cnt++;
    }
  }

  // Read the data from the file
  cnt = 0;
  //while (std::getline(iss, line)) {
  for (n=0;n<nscal;n++){
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
        for (i=0;i<nx;i++){
          std::getline(iss, line);
          std::istringstream sinput(line);
          sinput >> data[cnt];
          cnt++;
        }
      }
    }
  }
  //if(cnt != PeleC::h_prob_parm_device->nscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny){
  //  printf("Number of lines = %i. Specified in the header = %i \n", cnt,PeleC::h_prob_parm_device->nscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny);
  //  amrex::Abort("Number of lines in file differ from specified in header");
  //}
};

void
get_inputs(
  const std::string& iname,
  int z_coord)
{
  std::ifstream infile(iname, std::ios::in);
  const std::string memfile = read_file(infile); //The problem is coming from here
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }
  infile.close();
  std::istringstream iss(memfile);

  // Read the file
  std::string line;
  
  std::getline(iss, line);
  std::istringstream ninput(line);
  ninput >> PeleC::prob_parm_host->nscal_tmp;

  std::getline(iss, line);
  std::istringstream xinput(line);
  xinput >> PeleC::prob_parm_host->nx_tmp;
  
  std::getline(iss, line);
  std::istringstream yinput(line);
  yinput >> PeleC::prob_parm_host->ny_tmp;

  amrex::Print() << "Reading data from: " <<  iname.c_str() << std::endl;
  if(z_coord){
    std::getline(iss, line);
    std::istringstream zinput(line);
    zinput >> PeleC::prob_parm_host->nz_tmp;
    PeleC::prob_parm_host->nscal_tmp = PeleC::prob_parm_host->nscal_tmp - 3; //remove x and y and z from count
    amrex::Print() << PeleC::prob_parm_host->nscal_tmp << " scalars found with nx =  " << PeleC::prob_parm_host->nx_tmp 
                          << " and ny = " << PeleC::prob_parm_host->ny_tmp << " and nz = " << PeleC::prob_parm_host->nz_tmp << " points \n" <<std::endl;
  }
  else{
    PeleC::prob_parm_host->nscal_tmp = PeleC::prob_parm_host->nscal_tmp - 2; //remove x and y from count
    PeleC::prob_parm_host->nz_tmp = 1;
    amrex::Print() << PeleC::prob_parm_host->nscal_tmp << " scalars found with nx =  " << PeleC::prob_parm_host->nx_tmp 
                          << " and ny = " << PeleC::prob_parm_host->ny_tmp << " points \n" << std::endl;
  }

};


void PeleC::problem_post_timestep()
{
}

void PeleC::problem_post_init(//){}
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata)
{
}

void PeleC::problem_post_restart(){}
//  int i,
//  int j,
//  int k,
//  amrex::Array4<amrex::Real> const& state,
//  amrex::GeometryData const& geomdata){ 
//  int i,
//  int j,
//  int k,
//  amrex::Array4<amrex::Real> const& state,
//  amrex::GeometryData const& geomdata)
//{

 //  amrex::Real P,e;
 //  amrex::Real T = state(i,j,k,UTEMP);
 //  amrex::Real rho = state(i,j,k,URHO);
 //  amrex::Real massfrac[NUM_SPECIES];
 //  amrex::Real velx, vely, velz;
 //  for (int n = 0; n < NUM_SPECIES; n++)
 //        massfrac[n] = state(i, j, k, UFS + n)/rho;

 //  auto eos = pele::physics::PhysicsType::eos();
 //  eos.RTY2P(rho, T, massfrac, P);
 //  if(P > 4.5e+6){
 //    printf("Initial quantities:\n");
 //    printf("P[bar] = %f \n",P/1.0e+6);
 //    printf("T = %f \n",T);
 //    printf("U = %f \n",state(i, j, k, UMX)/state(i, j, k, URHO));
 //    printf("V = %f \n",state(i, j, k, UMY)/state(i, j, k, URHO));
 //    printf("W = %f \n",state(i, j, k, UMZ)/state(i, j, k, URHO));
 //    velx = 0.2*state(i, j, k, UMX)/rho;
 //    vely = 0.2*state(i, j, k, UMY)/rho;
 //    velz = 0.2*state(i, j, k, UMZ)/rho;
 //   //  //take the average velocity from neighboring cells
 //   //  if(i>0 and j>0 and k>0){

 //   //    velx = (state(i+1, j, k, UMX)/state(i+1, j, k, URHO) + state(i-1, j, k, UMX)/state(i-1, j, k, URHO) +
 //   //                      state(i, j+1, k, UMX)/state(i, j+1, k, URHO) + state(i, j-1, k, UMX)/state(i, j-1, k, URHO))/4.0;

 //   //    vely = (state(i+1, j, k, UMY)/state(i+1, j, k, URHO) + state(i-1, j, k, UMY)/state(i-1, j, k, URHO) +
 //   //                      state(i, j+1, k, UMY)/state(i, j+1, k, URHO) + state(i, j-1, k, UMY)/state(i, j-1, k, URHO))/4.0;

 //   //    velz = (state(i+1, j, k, UMZ)/state(i+1, j, k, URHO) + state(i-1, j, k, UMZ)/state(i-1, j, k, URHO) +
 //   //                      state(i, j+1, k, UMZ)/state(i, j+1, k, URHO) + state(i, j-1, k, UMZ)/state(i, j-1, k, URHO))/4.0;

 //   // //   velx = (state(i+1, j, k, UMX)/state(i+1, j, k, URHO) + state(i-1, j, k, UMX)/state(i-1, j, k, URHO) +
 //   // //                     state(i, j+1, k, UMX)/state(i, j+1, k, URHO) + state(i, j-1, k, UMX)/state(i, j-1, k, URHO) + 
 //   // //                     state(i, j, k+1, UMX)/state(i, j, k+1, URHO) + state(i, j, k-1, UMX)/state(i, j, k-1, URHO))/6.0;

 //   // //   vely = (state(i+1, j, k, UMY)/state(i+1, j, k, URHO) + state(i-1, j, k, UMY)/state(i-1, j, k, URHO) +
 //   // //                     state(i, j+1, k, UMY)/state(i, j+1, k, URHO) + state(i, j-1, k, UMY)/state(i, j-1, k, URHO) +   
 //   // //                     state(i, j, k+1, UMY)/state(i, j, k+1, URHO) + state(i, j, k-1, UMY)/state(i, j, k-1, URHO))/6.0;

 //   // //   velz = (state(i+1, j, k, UMZ)/state(i+1, j, k, URHO) + state(i-1, j, k, UMZ)/state(i-1, j, k, URHO) +
 //   // //                     state(i, j+1, k, UMZ)/state(i, j+1, k, URHO) + state(i, j-1, k, UMZ)/state(i, j-1, k, URHO) +   
 //   // //                     state(i, j, k+1, UMZ)/state(i, j, k+1, URHO) + state(i, j, k-1, UMZ)/state(i, j, k-1, URHO))/6.0;
 //   //    if(velx != velx or std::isinf(velx)){
 //   //       printf("velx = %f in i>0 and j>0 and k>0 i %i j %i k %i \n",velx,i,j,k);
 //   //       printf("state(i+1, j, k, UMX) = %f \n",state(i+1, j, k, URHO));
 //   //       printf("state(i-1, j, k, UMX) = %f \n",state(i-1, j, k, URHO));
 //   //       printf("state(i, j+1, k, UMX) = %f \n",state(i, j+1, k, URHO));
 //   //       printf("state(i, j-1, k, UMX) = %f \n",state(i, j-1, k, URHO));
 //   //    }
 //   //    if(vely != vely or std::isinf(vely)) printf("vely = %f in i>0 and j>0 and k>0 i %i j %i k %i \n",velz,i,j,k);
 //   //    if(velz != velz or std::isinf(velz)) printf("velz = %f in i>0 and j>0 and k>0 i %i j %i k %i \n",velz,i,j,k);

 //   //  }
 //   //  else if(i==0 and j>0 and k>0){
 //   //    velx = (state(i+1, j, k, UMX)/state(i+1, j, k, URHO) +
 //   //                      state(i, j+1, k, UMX)/state(i, j+1, k, URHO) + state(i, j-1, k, UMX)/state(i, j-1, k, URHO) +
 //   //                      state(i, j, k+1, UMX)/state(i, j, k+1, URHO) + state(i, j, k-1, UMX)/state(i, j, k-1, URHO))/5.0;

 //   //    vely = (state(i+1, j, k, UMY)/state(i+1, j, k, URHO) +
 //   //                      state(i, j+1, k, UMY)/state(i, j+1, k, URHO) + state(i, j-1, k, UMY)/state(i, j-1, k, URHO) +
 //   //                      state(i, j, k+1, UMY)/state(i, j, k+1, URHO) + state(i, j, k-1, UMY)/state(i, j, k-1, URHO))/5.0;

 //   //    velz = (state(i+1, j, k, UMZ)/state(i+1, j, k, URHO) +
 //   //                      state(i, j+1, k, UMZ)/state(i, j+1, k, URHO) + state(i, j-1, k, UMZ)/state(i, j-1, k, URHO) +
 //   //                      state(i, j, k+1, UMZ)/state(i, j, k+1, URHO) + state(i, j, k-1, UMZ)/state(i, j, k-1, URHO))/5.0;
 //   //    if(velx != velx or std::isinf(velx)) printf("velx = %f in i==0 and j>0 and k>0 \n",velx);
 //   //    if(vely != vely or std::isinf(vely)) printf("vely = %f in i==0 and j>0 and k>0 \n",velz);
 //   //    if(velz != velz or std::isinf(velz)) printf("velz = %f in i==0 and j>0 and k>0 \n",velz);
 //   //  }
 //   //  else if(i>0 and j==0 and k>0){
 //   //    velx = (state(i+1, j, k, UMX)/state(i+1, j, k, URHO) + state(i-1, j, k, UMX)/state(i-1, j, k, URHO) +
 //   //                      state(i, j+1, k, UMX)/state(i, j+1, k, URHO) +
 //   //                      state(i, j, k+1, UMX)/state(i, j, k+1, URHO) + state(i, j, k-1, UMX)/state(i, j, k-1, URHO))/5.0;

 //   //    vely = (state(i+1, j, k, UMY)/state(i+1, j, k, URHO) + state(i-1, j, k, UMY)/state(i-1, j, k, URHO) +
 //   //                      state(i, j+1, k, UMY)/state(i, j+1, k, URHO) +
 //   //                      state(i, j, k+1, UMY)/state(i, j, k+1, URHO) + state(i, j, k-1, UMY)/state(i, j, k-1, URHO))/5.0;

 //   //    velz = (state(i+1, j, k, UMZ)/state(i+1, j, k, URHO) + state(i-1, j, k, UMZ)/state(i-1, j, k, URHO) +
 //   //                      state(i, j+1, k, UMZ)/state(i, j+1, k, URHO) +
 //   //                      state(i, j, k+1, UMZ)/state(i, j, k+1, URHO) + state(i, j, k-1, UMZ)/state(i, j, k-1, URHO))/5.0;
 //   //    if(velx != velx or std::isinf(velx)) printf("velx = %f in i>0 and j==0 and k>0 \n",velx);
 //   //    if(vely != vely or std::isinf(vely)) printf("vely = %f in i>0 and j==0 and k>0 \n",velz);
 //   //    if(velz != velz or std::isinf(velz)) printf("velz = %f in i>0 and j==0 and k>0 \n",velz);
 //   //  }
 //   //  else if(i>0 and j>0 and k==0){

 //   //    velx = (state(i+1, j, k, UMX)/state(i+1, j, k, URHO) + state(i-1, j, k, UMX)/state(i-1, j, k, URHO) +
 //   //                      state(i, j+1, k, UMX)/state(i, j+1, k, URHO) + state(i, j-1, k, UMX)/state(i, j-1, k, URHO) +
 //   //                      state(i, j, k+1, UMX)/state(i, j, k+1, URHO) )/5.0;

 //   //    vely = (state(i+1, j, k, UMY)/state(i+1, j, k, URHO) + state(i-1, j, k, UMY)/state(i-1, j, k, URHO) +
 //   //                      state(i, j+1, k, UMY)/state(i, j+1, k, URHO) + state(i, j-1, k, UMY)/state(i, j-1, k, URHO) +
 //   //                      state(i, j, k+1, UMY)/state(i, j, k+1, URHO) )/5.0;

 //   //    velz = (state(i+1, j, k, UMZ)/state(i+1, j, k, URHO) + state(i-1, j, k, UMZ)/state(i-1, j, k, URHO) +
 //   //                      state(i, j+1, k, UMZ)/state(i, j+1, k, URHO) + state(i, j-1, k, UMZ)/state(i, j-1, k, URHO) +
 //   //                      state(i, j, k+1, UMZ)/state(i, j, k+1, URHO) )/5.0;
 //   //    if(velx != velx or std::isinf(velx)) printf("velx = %f in i>0 and j>0 and k==0 \n",velx);
 //   //    if(vely != vely or std::isinf(vely)) printf("vely = %f in i>0 and j>0 and k==0 \n",velz);
 //   //    if(velz != velz or std::isinf(velz)) printf("velz = %f in i>0 and j>0 and k==0 \n",velz);
 //   //  }

 //    P = 4.5e+6;
 //    if(T > T*0.3)T = T * 0.7;
 //    for (int n = 0; n < NUM_SPECIES; n++)
 //          massfrac[n] = state(i, j, k, UFS + n)/rho;

 //    eos.PYT2RE(P, massfrac, T, rho, e);

 //    state(i,j,k,URHO) = rho;
 //    state(i,j,k,UTEMP) = T;
 //    state(i, j, k, UMX) = velx * rho;
 //    state(i, j, k, UMY) = vely * rho;
 //    state(i, j, k, UMZ) = velz * rho;
 //    state(i, j, k, UEINT) = rho * e;
 //    state(i, j, k, UEDEN) = rho*(e+0.5*(velx*velx + vely*vely + velz*velz));
 // //   state(i, j, k, UTEMP) = T;
 //    for (int n = 0; n < NUM_SPECIES; n++)
 //        state(i, j, k, UFS + n) = massfrac[n] * rho;
 //  printf("Final state:\n");
 //  printf("Density = %f \n",state(i,j,k,URHO));
 //  printf("velx    = %f \n",state(i,j,k,UMX)/rho);
 //  printf("vely    = %f \n",state(i,j,k,UMY)/rho);
 //  printf("velz    = %f \n",state(i,j,k,UMZ)/rho);
 //  printf("internal energy = %f \n",state(i, j, k, UEINT));
 //  printf("total energy = %f \n",state(i, j, k, UEDEN));
 //  }


// void PeleC::problem_post_restart(  
//   int i,
//   int j,
//   int k,
//   amrex::Array4<amrex::Real> const& state,
//   amrex::GeometryData const& geomdata)
// {
  // const amrex::Real* prob_lo = geomdata.ProbLo();
  // const amrex::Real* prob_hi = geomdata.ProbHi();
  // const amrex::Real* dx = geomdata.CellSize();
  // const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  // const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  // const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  // amrex::Real u[3] = {0.0};
  // amrex::Real molefrac[NUM_SPECIES] = {0.0};
  // amrex::Real massfrac[NUM_SPECIES] = {0.0};
  // amrex::Real rho, energy, T, pres, e;

  // amrex::Real kernel_location[3];
  // const amrex::Real kernel_diameter = 100.0e-04;
  // const amrex::Real kernel_y_loc    = 1.92;
  // const amrex::Real kernel_height   = 1.01e-01;
  // double r;

  // kernel_location[0] = 0.0;
  // kernel_location[1] = kernel_y_loc-kernel_height/2.0;
  // kernel_location[2] = 0.0;

  // if(PeleC::h_prob_parm_device->init_kernel){
    
  //   r = sqrt(pow((kernel_location[0]-x),2)+pow((kernel_location[2]-z),2));
  //   if(y > h_prob_parm_device->yarray[0] and y < h_prob_parm_device->yarray[PeleC::h_prob_parm_device->ny-1] 
  //     and r < PeleC::h_prob_parm_device->xarray[PeleC::h_prob_parm_device->nx-1]){

  //     // Fill in the velocities and interpolated quantities.
  //     amrex::Real u[3] = {0.0};
  //     amrex::Real uinterp[PeleC::h_prob_parm_device->nscal];

  //     // Interpolation factors
  //     amrex::Real mod[3] = {0.0};
  //     int idx   = 0;
  //     int idxp1 = 0;
  //     int idy   = 0;
  //     int idyp1 = 0;
  //     amrex::Real slp[3] = {0.0};
      
  //     mod[0] = std::fmod(r, PeleC::h_prob_parm_device->Lxinput);
  //     locate(PeleC::h_prob_parm_device->xarray, PeleC::h_prob_parm_device->nx, mod[0], idx);
  //     idxp1 = (idx + 1) % PeleC::h_prob_parm_device->nx;
  //     //idxp1 = (idx + 1);
  //     slp[0] = (mod[0] - PeleC::h_prob_parm_device->xarray[idx]) / PeleC::h_prob_parm_device->xdiff[idx];

  //     mod[1] = std::fmod(y, PeleC::h_prob_parm_device->Lyinput);
  //     locate(h_prob_parm_device->yarray, PeleC::h_prob_parm_device->ny, mod[1], idy);
  //     idyp1 = (idy + 1) % PeleC::h_prob_parm_device->ny;
  //     //idyp1 = (idy + 1);
  //     slp[1] = (mod[1] - h_prob_parm_device->yarray[idy]) / PeleC::h_prob_parm_device->ydiff[idy];

  //     const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]);
  //     const amrex::Real f1 = slp[0] * (1 - slp[1]);
  //     const amrex::Real f2 = (1 - slp[0]) * slp[1];
  //     const amrex::Real f3 = slp[0] * slp[1];

  //     // Interpolate data
  //     for (int iscal=0;iscal<PeleC::h_prob_parm_device->nscal;iscal++){

  //       uinterp[iscal] = PeleC::h_prob_parm_device->data_ic[idx  +idy  *PeleC::h_prob_parm_device->nx+iscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny]*f0 +
  //                        PeleC::h_prob_parm_device->data_ic[idxp1+idy  *PeleC::h_prob_parm_device->nx+iscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny]*f1 +
  //                        PeleC::h_prob_parm_device->data_ic[idx  +idyp1*PeleC::h_prob_parm_device->nx+iscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny]*f2 +
  //                        PeleC::h_prob_parm_device->data_ic[idxp1+idyp1*PeleC::h_prob_parm_device->nx+iscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny]*f3;
      
  //     }

  //     // // temperature, pressure, internal energy and composition
  //     double p ;//= uinterp[1]*10.0; //Pa to dyn/cm2
  //     rho = uinterp[0]/1.0e+3;
  //     T = uinterp[2];
  //     double massfrac[NUM_SPECIES]={0.0};

  //     for (int n = 0; n < NUM_SPECIES; n++)
  //       massfrac[n] = uinterp[n+6];

  //     double sum = 0.0;
  //     for (int n = 0; n < NUM_SPECIES; n++)
  //       sum = sum + massfrac[n];
  //     // printf("Sum of all species = % e \n", sum);
  //     massfrac[N2_ID] = massfrac[N2_ID] + 1.0 - sum; //making sure that sum of ys is 1


  //     auto eos = pele::physics::PhysicsType::eos();
  //     eos.RTY2P(rho, T, massfrac, p);
  //     eos.PYT2RE(p, massfrac, T, rho, e);

  //     double beta  = atan(z/x);
  //     double alpha = atan(x/z);

  //     //Velocity in each direction
  //     if(x >= 0.0){
  //       u[0] = uinterp[5]*cos(beta)*100.0;
  //     }
  //     else{
  //       u[0] = -uinterp[5]*cos(beta)*100.0;
  //     }

  //     u[1] = uinterp[4]*100.0;

  //     if(z >= 0){
  //       u[2] = uinterp[5]*cos(alpha)*100.0;
  //     }
  //     else{
  //       u[2] = -uinterp[5]*cos(alpha)*100.0;
  //     }
  //     // for (int n = 0; n < NUM_SPECIES; n++)
  //     //   state(i, j, k, UFS + n) = rho * massfrac[n];

  //     amrex::Real rho_init = state(i, j, k, URHO);
  //     amrex::Real velx = uinterp[0];
  //     amrex::Real vely = uinterp[1];
  //     amrex::Real velz = uinterp[2];
  //     // printf("rho_init = %f \n", rho_init);
  //     // printf("UMX      = %f\n", state(i, j, k, UMX)/rho_init);
  //     // printf("UMY      = %f\n", state(i, j, k, UMY)/rho_init);
  //     // printf("UMZ      = %f\n", state(i, j, k, UMZ)/rho_init);

  //     // for (int n = 0; n < NUM_SPECIES; n++)
  //     //   state(i, j, k, UFS + n) = rho_init * massfrac[n];

  //     // state(i, j, k, URHO)  = rho;
  //     // state(i, j, k, UEINT) = rho * e;
  //     // state(i, j, k, UTEMP) = T;
  //     // state(i, j, k, UMX)   = rho * velx;
  //     // state(i, j, k, UMY)   = rho * vely;
  //     // state(i, j, k, UMZ)   = rho * velz;
  //     // state(i, j, k, UEDEN)  = rho * (e + 0.5 * (velx * velx + vely * vely + velz * velz)); 

  //     if(state(i, j, k, URHO) != state(i, j, k, URHO))
  //       printf("state(i, j, k, URHO) = %f\n", state(i, j, k, URHO)); 

  //     if(state(i, j, k, UEINT) != state(i, j, k, UEINT))
  //       printf("state(i, j, k, UEINT) = %f\n", state(i, j, k, UEINT)); 

  //     if(state(i, j, k, UTEMP) != state(i, j, k, UTEMP))
  //       printf("state(i, j, k, UTEMP) = %f\n", state(i, j, k, UTEMP)); 

  //     if(state(i, j, k, UMX) != state(i, j, k, UMX))
  //       printf("state(i, j, k, UMX) = %f\n", state(i, j, k, UMX)); 

  //     if(state(i, j, k, UMY) != state(i, j, k, UMY))
  //       printf("state(i, j, k, UMY) = %f\n", state(i, j, k, UMY)); 

  //     if(state(i, j, k, UMZ) != state(i, j, k, UMZ))
  //       printf("state(i, j, k, UMZ) = %f\n", state(i, j, k, UMZ)); 

  //     if(state(i, j, k, UEDEN) != state(i, j, k, UEDEN))
  //       printf("state(i, j, k, UEDEN) = %f\n", state(i, j, k, UEDEN)); 

  //   }      
  //   // else{
  //   //   for (int n = 0; n < NVAR; n++)
  //   //     state(i, j, k, n) = PeleC::h_prob_parm_device->fuel_state[n];
  //   // }
  // }

//}
