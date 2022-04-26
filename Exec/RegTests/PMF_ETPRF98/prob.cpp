#include "prob.H"

std::string
read_pmf_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

bool
checkQuotes(const std::string& str)
{
  int count = 0;
  for (char c : str) {
    if (c == '"') {
      count++;
    }
  }
  return (count % 2) == 0;
}

void
read_pmf(const std::string& myfile)
{
  std::string firstline;
  std::string secondline;
  std::string remaininglines;
  unsigned int pos1;
  unsigned int pos2;
  int variable_count;
  int line_count;

  std::ifstream infile(myfile);
  const std::string memfile = read_pmf_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  if (!checkQuotes(firstline)) {
    amrex::Abort("PMF file variable quotes unbalanced");
  }
  std::getline(iss, secondline);
  pos1 = 0;
  pos2 = 0;
  variable_count = 0;
  while ((pos1 < firstline.length() - 1) && (pos2 < firstline.length() - 1)) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    variable_count++;
    pos1 = pos2 + 1;
  }

  amrex::Vector<std::string> pmf_names;
  pmf_names.resize(variable_count);
  pos1 = 0;
  // pos2 = 0;
  for (int i = 0; i < variable_count; i++) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    pmf_names[i] = firstline.substr(pos1 + 1, pos2 - (pos1 + 1));
    pos1 = pos2 + 1;
  }

  amrex::Print() << variable_count << " variables found in PMF file"
                 << std::endl;
  // for (int i = 0; i < variable_count; i++)
  //  amrex::Print() << "Variable found: " << pmf_names[i] <<
  //  std::endl;

  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in PMF file" << std::endl;

  PeleC::h_prob_parm_device->pmf_N = line_count;
  PeleC::h_prob_parm_device->pmf_M = variable_count - 1;
  PeleC::prob_parm_host->h_pmf_X.resize(PeleC::h_prob_parm_device->pmf_N);
  PeleC::prob_parm_host->pmf_X.resize(PeleC::h_prob_parm_device->pmf_N);
  PeleC::prob_parm_host->h_pmf_Y.resize(
    PeleC::h_prob_parm_device->pmf_N * PeleC::h_prob_parm_device->pmf_M);
  PeleC::prob_parm_host->pmf_Y.resize(
    PeleC::h_prob_parm_device->pmf_N * PeleC::h_prob_parm_device->pmf_M);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  std::getline(iss, secondline);
  for (unsigned int i = 0; i < PeleC::h_prob_parm_device->pmf_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> PeleC::prob_parm_host->h_pmf_X[i];
    for (unsigned int j = 0; j < PeleC::h_prob_parm_device->pmf_M; j++) {
      sinput >> PeleC::prob_parm_host
                  ->h_pmf_Y[j * PeleC::h_prob_parm_device->pmf_N + i];
    }
  }

  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_X.begin(),
    PeleC::prob_parm_host->h_pmf_X.end(), PeleC::prob_parm_host->pmf_X.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_Y.begin(),
    PeleC::prob_parm_host->h_pmf_Y.end(), PeleC::prob_parm_host->pmf_Y.begin());
  PeleC::h_prob_parm_device->d_pmf_X = PeleC::prob_parm_host->pmf_X.data();
  PeleC::h_prob_parm_device->d_pmf_Y = PeleC::prob_parm_host->pmf_Y.data();
  PeleC::d_prob_parm_device->d_pmf_X = PeleC::prob_parm_host->pmf_X.data();
  PeleC::d_prob_parm_device->d_pmf_Y = PeleC::prob_parm_host->pmf_Y.data();
}

void
init_bc()
{
  amrex::Real vt;
  amrex::Real ek;
  amrex::Real T;
  amrex::Real rho;
  amrex::Real e;
  amrex::Real molefrac[NUM_SPECIES];
  amrex::Real massfrac[NUM_SPECIES];
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};

  if (PeleC::h_prob_parm_device->phi_in < 0) {
    const amrex::Real yl = 0.0;
    const amrex::Real yr = 0.0;
    pmf(yl, yr, pmf_vals, *PeleC::h_prob_parm_device);
    amrex::Real mysum = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      molefrac[n] = amrex::max<amrex::Real>(0.0, pmf_vals[3 + n]);
      mysum += molefrac[n];
    }
    molefrac[N2_ID] = 1.0 - (mysum - molefrac[N2_ID]);
    T = pmf_vals[0];
    PeleC::h_prob_parm_device->vn_in = pmf_vals[1];
  } else {
    
    //number of carbon and hydrogen atoms in each fuel component
    //order: IC8H18; NC7H16_ID; C6H5CH3; C2H5OH
    amrex::Real C_atoms[4] = {8.0 ,7.0 ,7.0,2.0}; 
    amrex::Real H_atoms[4] = {18.0,16.0,4.0,6.0}; 
    amrex::Real O_atoms[4] = {0.0 ,0.0 ,0.0,1.0}; 

    amrex::Real volfrac_IC8H18  = 0.284;
    amrex::Real volfrac_NC7H16  = 0.172;
    amrex::Real volfrac_C6H5CH3 = 0.339;
    amrex::Real volfrac_C2H5OH  = 0.205;
     
    amrex::Real volfrac_CO2 =  C_atoms[0]*volfrac_IC8H18 + C_atoms[1]*volfrac_NC7H16 + C_atoms[2]*volfrac_C6H5CH3 + C_atoms[3]*volfrac_C2H5OH;
    amrex::Real volfrac_H2O = (H_atoms[0]*volfrac_IC8H18 + H_atoms[1]*volfrac_NC7H16 + H_atoms[2]*volfrac_C6H5CH3 + H_atoms[3]*volfrac_C2H5OH)/2.0;
     
    amrex::Real volfrac_O2  = (2.0*volfrac_CO2+volfrac_H2O-(O_atoms[0]*volfrac_IC8H18  + 
                                                            O_atoms[1]*volfrac_NC7H16  + 
                                                            O_atoms[2]*volfrac_C6H5CH3 + 
                                                            O_atoms[3]*volfrac_C2H5OH))/2.0/PeleC::h_prob_parm_device->phi_in;
    amrex::Real volfrac_N2  = volfrac_O2*(79./21.);

    amrex::Real sum_oxi  = volfrac_IC8H18 + volfrac_NC7H16 + volfrac_C6H5CH3 + volfrac_C2H5OH + volfrac_O2  + volfrac_N2;
    amrex::Real sum_prod = volfrac_CO2    + volfrac_H2O    + volfrac_N2;

    amrex::Real X_IC8H18  = volfrac_IC8H18   / sum_oxi;
    amrex::Real X_NC7H16  = volfrac_NC7H16   / sum_oxi;
    amrex::Real X_C6H5CH3 = volfrac_C6H5CH3  / sum_oxi;
    amrex::Real X_C2H5OH  = volfrac_C2H5OH   / sum_oxi;
    amrex::Real X_O2      = volfrac_O2       / sum_oxi;
    amrex::Real X_N2      = volfrac_N2       / sum_oxi;

    amrex::Real X_co2    = volfrac_CO2 / sum_prod;
    amrex::Real X_h2o    = volfrac_H2O / sum_prod;

    for (int n = 0; n < NUM_SPECIES; n++)
      molefrac[n] = 0.0;

    molefrac[IC8_ID]     = X_IC8H18  * (1.0-PeleC::h_prob_parm_device->egr);
    molefrac[NC7H16_ID]  = X_NC7H16  * (1.0-PeleC::h_prob_parm_device->egr);
    molefrac[C6H5CH3_ID] = X_C6H5CH3 * (1.0-PeleC::h_prob_parm_device->egr);
    molefrac[C2H5OH_ID]  = X_C2H5OH  * (1.0-PeleC::h_prob_parm_device->egr);
    molefrac[O2_ID]      = X_O2      * (1.0-PeleC::h_prob_parm_device->egr);
    molefrac[N2_ID]      = X_N2      * (1.0-PeleC::h_prob_parm_device->egr) + (PeleC::h_prob_parm_device->egr*volfrac_N2 / sum_prod);
    molefrac[CO2_ID]     = X_co2*PeleC::h_prob_parm_device->egr;
    molefrac[H2O_ID]     = X_h2o*PeleC::h_prob_parm_device->egr;

    T = PeleC::h_prob_parm_device->T_in;
  }
  const amrex::Real p = PeleC::h_prob_parm_device->pamb;

  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(molefrac, massfrac);
  eos.PYT2RE(p, massfrac, T, rho, e);

  vt = PeleC::h_prob_parm_device->vn_in;
  ek = 0.5 * (vt * vt);

  PeleC::h_prob_parm_device->fuel_state[URHO] = rho;
  PeleC::h_prob_parm_device->fuel_state[UMX] = 0.0;
  PeleC::h_prob_parm_device->fuel_state[UMY] = rho * vt;
  PeleC::h_prob_parm_device->fuel_state[UMZ] = 0.0;
  PeleC::h_prob_parm_device->fuel_state[UEINT] = rho * e;
  PeleC::h_prob_parm_device->fuel_state[UEDEN] = rho * (e + ek);
  PeleC::h_prob_parm_device->fuel_state[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    PeleC::h_prob_parm_device->fuel_state[UFS + n] = rho * massfrac[n];
  }
}

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  std::string pmf_datafile;

  amrex::ParmParse pp("prob");
  pp.query("pamb", PeleC::h_prob_parm_device->pamb);
  pp.query("phi_in", PeleC::h_prob_parm_device->phi_in);
  pp.query("T_in", PeleC::h_prob_parm_device->T_in);
  pp.query("vn_in", PeleC::h_prob_parm_device->vn_in);
  pp.query("pertmag", PeleC::h_prob_parm_device->pertmag);
  pp.query("pmf_datafile", pmf_datafile);

  PeleC::h_prob_parm_device->L[0] = probhi[0] - problo[0];
  PeleC::h_prob_parm_device->L[1] = probhi[1] - problo[1];
  PeleC::h_prob_parm_device->L[2] = probhi[2] - problo[2];

  read_pmf(pmf_datafile);

  init_bc();
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init(){}
//  int i,
//  int j,
//  int k,
//  amrex::Array4<amrex::Real> const& state,
//  amrex::GeometryData const& geomdata)
//{
//}

void
PeleC::problem_post_restart()
{
}
