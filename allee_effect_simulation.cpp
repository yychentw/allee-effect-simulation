#include <iostream>
#include <math.h>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
using namespace std;

//To generate random number
class RandGenerator{
public:
	RandGenerator();
	//to generate a random number from the continuous uniform distribution on the interval [start, end]
    double	RUnif(double Start, double End);
	//to select a random item from the list {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0} to determine an individual's degree of cooperation
	double	RandDeg();
	//to generate a random number from a Bernoulli distribution
    bool	RBern(double Prob);
	//to generate a random number from a 2D Gaussian distribution
	int		R2DGaussian();
	//to generate a random number from a uniform integer distribution on the interval [start, end-1]
	int		SampleInt(int Start, int End);

private:
	unsigned seed;
	default_random_engine RAND_GENERATOR;
};

class Individual{
public:
    Individual(int PatchID, double CoopDeg, int Age);
	void	Aging(){ _Age = _Age + 1; }
    void	UpdateReproRate(int    ReproRate){ _ReproRate = ReproRate; }
	void	UpdateSurvRate (double SurvRate ){ _SurvRate  = SurvRate;  }
	int		GetAge()        { return _Age;         }
    int		GetReproRate()  { return _ReproRate;   }
	double	GetSurvRate()   { return _SurvRate;    }
    double	GetCoopDeg()    { return _CoopDeg;     }

private:
	int		_Age;
    int     _PatchID;		// ID of the patch to which the individual belongs
    int     _ReproRate;		// reproduction rate
	double	_SurvRate;		// survival rate
    double  _CoopDeg;		// degree of cooperation
	RandGenerator RG;
};

class Patch{
public:
	Patch();
	// to add an individual to the patch
    void	AddIndividual(int PatchID, double CoopDeg, int Age);
    void	UpdateCoopNum();
    void	UpdateCoopDegSum();
	void	UpdateIndivAge();
    void	UpdateIndivReproRate(double CoopEfficiency, double HalfConst, double Temperature, double Topt, double CTmax, double ReproMax, double CostRate);
	void	UpdateIndivSurvRate (double SurvRateUL, double AgeStandard);
	void	UpdateIndivAlive();
	// to get the size of the patch (or the number of the individuals in the patch)
	int		GetSize()		   		  { return _Individual.size();  }
    int		GetCoopNum()   		      { return _CoopNum;            }
    double	GetCoopDegSum()		      { return _CoopDegSum;	        }
	// to get the total number of offspring will be produced by all the individuals in the patch
	vector<double> GetOffspring();

private:
    int     _CoopNum;		// the number of cooperators in the patch
    double  _CoopDegSum;	// the sum of the degree of cooperation of all individuals in the patch
	vector<Individual> _Individual;
	RandGenerator RG;
};

//To manage processes in population level
class PopProcess{
public:
	PopProcess(int OutputType, double HabitatAvail, double Temperature, double HabitatLength, double ResourceAvail, double InitPopDens, int Replic, int Span, \
			   double Topt, double CTmax, bool Sociality, double CoopEfficiency, double HalfConst, double ReproMax, double CostRate, double InitCoopProp, int RandomSeed);
	~PopProcess();
	// to compute the number of cooperators, the proportion of the cooperators, and average degree of cooperation at population level
    void	ComputePopCoopData();
	// to compute and print all information regarding the population status at the given time
    void    ComputeTimeSeriesData(int Time);
	// to handle cross-patches offspring dispersal
	void	OffspringDisperse();
	// to iteratively update the population until the EndTime
    void	SystemUpdate(int EndTime);
	// to handle the exception as the population go extinct
	void	ErrorHandling1();
	// to handle the exception as the population becomes larger than the upper limit
	void	ErrorHandling2();
	int		GetPopSize()				{ return _PopSize;	        }
	int		GetOffspringNum()			{ return _OffspringNum;		}
    int		GetTotalCoopNum()   		{ return _TotalCoopNum;	    }
    double	GetCoopProp()   			{ return _CoopProp;		    }
    double	GetAverCoopDeg()    		{ return _AverCoopDeg;	    }
	vector<bool>	GetPatchIsAvail()	{ return _PatchIsAvail;		}
	vector<Patch>	GetPatch()			{ return _Patch;			}

private:
	//Population state
    int     _PopSize;			// population size
	int		_PatchNum;			// the number of patches in the population
    int		_PopSizeUL;			// the upper limit of the population size
    int		_TotalCoopNum;		// the total number of cooperators in the population
	double	_InitPopDens;		// initial population density
    double	_CoopProp;			// the proportion of cooperators in the population
	double	_InitCoopProp;		// the initial proportion of cooperators in the population
    double	_AverCoopDeg;		// average degree of cooperation
	bool	_Sociality;			// the sociality of the population, 1 means a social population and 0 means a non-social population 
	
    //Population structure
	vector<Patch> _Patch;

	//Environmental factor
	double	_HabitatAvail;		// the availability of suitable habitat
	double	_ResourceAvail;		// the availability of resources
	double	_Temperature;		// environmental temperature

	//Reproduction
	double	_Topt;				// the optimal temperature for reproduction
	double	_CTmax;				// the critical maximum temperature for reproduction
    double	_CoopEfficiency;	// the efficiency of transfering cooperation efforts to cooperation benefits
    double	_CostRate;			// the percentage decrease in the reproduction rate caused by per unit cooperation degree
    double  _HalfConst;			// the “half-saturation constant”, which is the value of the cooperation benefits at which the reproduction gain is half of its maximum
    double  _ReproMax;			// the maximum reproduction rate without cooperation
	int		_OffspringNum;		// the number of offspring in the population

	//Survival rate
	double	_AgeStandard;		// the exponential age constant to determine the exponential decay of the survival rate
    double	_SurvRateUL;		// the upper limit of survival rate
	
    //Random number generator
	RandGenerator RG;

    //For convenience of data computation
    int		_Span;				// the time interval which determines the frequency of the time series output
	int		_Replic;			// the serial number of the simulation replication
	int		_Length;			// the length of the square habitat
	int		_OutputType;
	vector<int> _AvailPatchID;	// the IDs of the patches which are in the suitable habitat
	vector<bool>_PatchIsAvail;
};



void PopProcess::SystemUpdate(int EndTime){
	for (int time = 0; time < EndTime; time++){

		// to select patches which are with enough resources for breeding
		vector<int> breeding_patch_id;
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		default_random_engine generator(seed);
		shuffle(_AvailPatchID.begin(), _AvailPatchID.end(), generator);
		for (int i = 0; i < int(_AvailPatchID.size()*_ResourceAvail); i++){
			breeding_patch_id.push_back(_AvailPatchID[i]);
		}

        // to update all patches in the population
		for (Patch& patch: _Patch){
			patch.UpdateIndivAge();
            patch.UpdateCoopNum();
            patch.UpdateCoopDegSum();
		}
		for (int patch_id: breeding_patch_id){
			_Patch[patch_id].UpdateIndivReproRate(_CoopEfficiency, _HalfConst, _Temperature, _Topt, _CTmax, _ReproMax, _CostRate);
		}
		OffspringDisperse();
		for (int i = 0; i < _PatchNum; i++){
			_Patch[i].UpdateIndivSurvRate(_SurvRateUL, _AgeStandard);
			_Patch[i].UpdateIndivAlive();
		}
		
        // to compute and check population size
		_PopSize = 0;
		for (Patch& patch: _Patch){
			_PopSize = _PopSize + patch.GetSize();
		}
		if (_PopSize == 0){
			ErrorHandling1();
            break;
		}
		if (_PopSize > _PopSizeUL){
			ErrorHandling2();
			break;
		}

		// to compute the number of cooperators, the proportion of the cooperators, and average degree of cooperation at population level
        ComputePopCoopData();

        // to compute time series data
		if ((time%_Span == _Span-1) && (_OutputType == 2)){
			ComputeTimeSeriesData(time);
		}
	}
}


class Simulation{
public:
	Simulation(int OutputType, double HabitatAvail, double Temperature, int HabitatLength, double ResourceAvail, double InitPopDens, int EndTime, int InitReplic, \
			   int ReplicNum, int Span, double Topt, double CTmax, bool Sociality, double CoopEfficiency, double HalfConst, double ReproMax, double CostRate, \
			   double InitCoopProp, int RandomSeed);
    void	ProcessShunt(int OutputType);
	void	PopTimeSeries();
	void	PopFinalState();

private:
	int		_OutputType;
	int		_HabitatLength;
	int		_EndTime;
	int		_InitReplic;
	int		_ReplicNum;
	int		_Span;
	int		_RandomSeed;
	bool	_Sociality;
	double	_InitPopDens;
	double	_HabitatAvail;
	double	_ResourceAvail;
	double	_Temperature;
	double	_Topt;
	double	_CTmax;
	double	_CoopEfficiency;
	double	_HalfConst;
	double	_ReproMax;
	double	_CostRate;
	double	_InitCoopProp;
};



int main(int arc, char *argv[]){

	int     	output_type, habitat_length, end_time, init_replic, replic_num, span, random_seed; // note that negative random seed will not be used
    bool    	sociality;
    double  	habitat_avail, temperature, topt, ctmax, resource_avail, init_pop_dens, coop_efficiency, half_const, repro_max, cost_rate, init_coop_prop;
	
	cin >> output_type >> habitat_avail >> habitat_length >> resource_avail >> init_pop_dens >> end_time >> init_replic >> replic_num >> span >> topt >> ctmax \
		>> sociality >> coop_efficiency >> half_const >> repro_max >> cost_rate >> init_coop_prop >> random_seed;
	
	int patch_num = int(pow(habitat_length, 2));
	
	// to print headers
	if (output_type == 1){
		cout << "replic,habitat_ayvail,temperature,habitat_length,resource_avail,init_pop_dens,topt,ctmax,coop_efficiency,half_const,repro_max,cost_rate,socialit,pop_size,offspring_num,coop_num,noncoop_num,coop_prop,aver_coop_deg,init_coop_prop";
	}else if (output_type == 2){
		cout << "replic,time,habitat_avail,temperature,habitat_length,resource_avail,init_pop_dens,topt,ctmax,coop_efficiency,half_const,repro_max,cost_rate,sociality,pop_size,offspring_num,coop_num,noncoop_num,coop_prop,aver_coop_deg,init_coop_prop";
	}
	for (int var_type = 0; var_type < 3; var_type++){
		for (int patch_id = 0; patch_id < patch_num; patch_id++){
			cout << ",";
			switch (var_type)
			{
				case 0:{
					cout << "forest_" << patch_id;
				}break;
				case 1:{
					cout << "group_size_" << patch_id;
				}break;
				case 2:{
					cout << "offspring_num_" << patch_id;
				}break;
			}
		}
	}
	cout << "\n";

	while (cin >> temperature){
		Simulation Sim(output_type, habitat_avail, temperature, habitat_length, resource_avail, init_pop_dens, end_time, init_replic, replic_num, span, topt, ctmax, \
					   sociality, coop_efficiency, half_const, repro_max, cost_rate, init_coop_prop, random_seed);
		Sim.ProcessShunt(output_type);
	}
    
	
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

RandGenerator::RandGenerator(){
	seed = chrono::system_clock::now().time_since_epoch().count();
	RAND_GENERATOR.seed(seed);
}

double RandGenerator::RUnif(double Start, double End){
	uniform_real_distribution<double> distribution(Start, End);
	return distribution(RAND_GENERATOR);
}

double RandGenerator::RandDeg(){
    uniform_int_distribution<int> distribution(1, 10);
    double output = (double)distribution(RAND_GENERATOR)/10;
    return output;
}

bool RandGenerator::RBern(double Prob){
	bernoulli_distribution distribution(Prob);
	return distribution(RAND_GENERATOR);
}

int RandGenerator::R2DGaussian(){
	discrete_distribution<int>  distribution {0.003765, 0.015019, 0.023792, 0.015019, 0.003765, \
                                           	  0.015019, 0.059912, 0.094907, 0.059912, 0.015019, \
                                           	  0.023792, 0.094907, 0.150342, 0.094907, 0.023792, \
                                           	  0.015019, 0.059912, 0.094907, 0.059912, 0.015019, \
                                           	  0.003765, 0.015019, 0.023792, 0.015019, 0.003765};
	return distribution(RAND_GENERATOR);
}

int RandGenerator::SampleInt(int Start, int End){
	uniform_int_distribution<int> distribution(Start, End-1);
	return distribution(RAND_GENERATOR);
}

Individual::Individual(int PatchID, double CoopDeg, int Age){
	_Age            = Age;
	_PatchID        = PatchID;
	_CoopDeg        = CoopDeg;
	_ReproRate      = 0;

	if (_CoopDeg == -1.0){
		//_CoopDeg = RG.RandDeg();
		_CoopDeg = 0.5;
	}
}

Patch::Patch(){
}

void Patch::AddIndividual(int PatchID, double CoopDeg, int Age){
	_Individual.push_back(Individual(PatchID, CoopDeg, Age));
}

void Patch::UpdateCoopDegSum(){
	_CoopDegSum = 0.0;
	for (Individual& indiv: _Individual){
		_CoopDegSum = _CoopDegSum + indiv.GetCoopDeg();
	}
}

void Patch::UpdateCoopNum(){
	_CoopNum = 0;
	for (Individual& indiv: _Individual){
		if (indiv.GetCoopDeg() > 0.0){
			_CoopNum++;
    	}
	}
}

void Patch::UpdateIndivAge(){
	for (Individual& indiv: _Individual){
		indiv.Aging();
	}
}

void Patch::UpdateIndivReproRate(double CoopEfficiency, double HalfConst, double Temperature, double Topt, double CTmax, double ReproMax, double CostRate){
	// to compute the effect of cooperation
	double coop_benefit = _CoopDegSum * CoopEfficiency * (Temperature - Topt)/(CTmax - Topt);
	double tmp = 1 + coop_benefit/(HalfConst + coop_benefit);

	double performance = 0;
	if (Temperature < CTmax){
		performance = 1 - pow((Temperature-Topt)/(Topt-CTmax), 2);
	}
	// to compute individual reproduction rate and total reproduction of the group (or patch)
	for (Individual& indiv: _Individual){
		double	expected_repro_rate;
		if (_Individual.size() > 1){
			expected_repro_rate = ReproMax/sqrt(_Individual.size()) * tmp * (1.0 - CostRate*indiv.GetCoopDeg()) * performance;
		}else{
			expected_repro_rate = ReproMax/sqrt(_Individual.size()) * performance;
		}

		if (expected_repro_rate <= 0){
			indiv.UpdateReproRate(0);
		}else{
			double  repro_rate_floor = floor(expected_repro_rate);
			double  prob = expected_repro_rate - repro_rate_floor;
			int     indiv_repro;
			if (RG.RBern(prob)){
				indiv_repro = (int)(repro_rate_floor+1);
			}else{
				indiv_repro = (int)repro_rate_floor;
			}
			indiv.UpdateReproRate(indiv_repro);
		}
	}
}

void Patch::UpdateIndivSurvRate(double SurvRateUL, double AgeStandard){
	for (Individual& indiv: _Individual){
		double IndivSurvRate = SurvRateUL*exp(-indiv.GetAge()/AgeStandard);
		indiv.UpdateSurvRate(IndivSurvRate);
	}
}

void Patch::UpdateIndivAlive(){
	for (int i = 0; i < _Individual.size(); i++){
		if (!RG.RBern(_Individual[i].GetSurvRate())){
			_Individual.erase(_Individual.begin()+i);
			i--;
		}
	}
}

vector<double> Patch::GetOffspring(){
	vector<double> Offspring;
	for (Individual& indiv: _Individual){
		vector<double> tmp(indiv.GetReproRate(), indiv.GetCoopDeg());
		Offspring.insert(Offspring.end(), tmp.begin(), tmp.end());
	}
	return Offspring;
}

PopProcess::PopProcess(int OutputType, double HabitatAvail, double Temperature, double HabitatLength, double ResourceAvail, double InitPopDens, int Replic, int Span, \
					   double Topt, double CTmax, bool Sociality, double CoopEfficiency, double HalfConst, double ReproMax, double CostRate, double InitCoopProp, 
					   int RandomSeed){
    //Population state
    _PatchNum			= int(pow(HabitatLength, 2));
    _PopSizeUL			= 1000000;
	_PopSize			= int(InitPopDens*_PatchNum);
	_InitPopDens		= InitPopDens;
	_InitCoopProp		= InitCoopProp;
	if (Sociality){
		_CoopProp		= InitCoopProp;
	}else{
		_CoopProp		= 0.0;
	}
	_Sociality			= Sociality;

    //Population structure
	_Patch = vector<Patch>(_PatchNum);

	//Environmental factor
	_HabitatAvail		= HabitatAvail;
	_ResourceAvail		= ResourceAvail;
	_Temperature		= Temperature;

	//Reproduction
	_Topt				= Topt;
	_CTmax				= CTmax;
    _CoopEfficiency		= CoopEfficiency;
    _CostRate			= CostRate;
    _HalfConst			= HalfConst;
    _ReproMax			= ReproMax;

	//Survival rate
	_AgeStandard		= 2.0;
    _SurvRateUL			= 0.7;

    //For convenience of data computation
	_Length				= HabitatLength;
    _Span				= Span;
	_Replic				= Replic;
	_OutputType			= OutputType;

	// to select the patches in the suitable habitat
	vector<int> tmp(_PatchNum, 0);
	for (int i = 0; i < _PatchNum; i++){
		tmp[i] = i;
	}
	unsigned seed;
	if (RandomSeed < 0){
		seed = chrono::system_clock::now().time_since_epoch().count();
	}else{
		seed = unsigned(RandomSeed);
	}
	default_random_engine generator(seed);
	shuffle(tmp.begin(), tmp.end(), generator);
	for (int i = 0; i < int(_PatchNum*_HabitatAvail); i++){
		_AvailPatchID.push_back(tmp[i]);
	}
	_PatchIsAvail = vector<bool>(_PatchNum, 0);
	for (int patch_id: _AvailPatchID){
		_PatchIsAvail[patch_id] = 1;
	}

	// to initialize population
	int patch_id;
	default_random_engine RAND_GENERATOR;
	RAND_GENERATOR.seed(seed);
	uniform_int_distribution<int> distribution(0, _PatchNum-1);
	for (int i = 0; i < _PopSize; i++){
		patch_id = distribution(RAND_GENERATOR);
		if (i < _PopSize*_CoopProp){
			_Patch[patch_id].AddIndividual(patch_id, -1.0, 0);
		}else{
			_Patch[patch_id].AddIndividual(patch_id, 0.0, 0);
		}
	}
	ComputePopCoopData();
}

PopProcess::~PopProcess(){}

void PopProcess::OffspringDisperse(){
	int patch_id;
	_OffspringNum = 0;
	for (int i = 0; i < _PatchNum; i++){
		vector<double> Offspring = _Patch[i].GetOffspring();
		_OffspringNum = _OffspringNum + Offspring.size();
		for (double parent_trait: Offspring){
			patch_id = i;
			int dir = RG.R2DGaussian();
			switch (dir){
				case 0:{
					if ((i >= _Length*2) && (i%_Length >= 2)){ patch_id = (i/_Length-2)*_Length + i%_Length-2;	}
				}break;
				case 1:{
					if ((i >= _Length*2) && (i%_Length >= 1)){ patch_id = (i/_Length-2)*_Length + i%_Length-1;	}
				}break;
				case 2:{
					if (i >= _Length*2){ patch_id = (i/_Length-2)*_Length + i%_Length;	}
				}break;
				case 3:{
					if ((i >= _Length*2) && (i%_Length <= _Length-2)){ patch_id = (i/_Length-2)*_Length + i%_Length+1;	}
				}break;
				case 4:{
					if ((i >= _Length*2) && (i%_Length <= _Length-3)){ patch_id = (i/_Length-2)*_Length + i%_Length+2;	}
				}break;
				case 5:{
					if ((i >= _Length) && (i%_Length >= 2)){ patch_id = (i/_Length-1)*_Length + i%_Length-2;	}
				}break;
				case 6:{
					if ((i >= _Length) && (i%_Length >= 1)){ patch_id = (i/_Length-1)*_Length + i%_Length-1;	}
				}break;
				case 7:{
					if (i >= _Length){ patch_id = (i/_Length-1)*_Length + i%_Length;	}
				}break;
				case 8:{
					if ((i >= _Length) && (i%_Length <= _Length-2)){ patch_id = (i/_Length-1)*_Length + i%_Length+1;	}
				}break;
				case 9:{
					if ((i >= _Length) && (i%_Length <= _Length-3)){ patch_id = (i/_Length-1)*_Length + i%_Length+2;	}
				}break;
				case 10:{
					if (i%_Length >= 2){ patch_id = i-2;	}
				}break;
				case 11:{
					if (i%_Length >= 1){ patch_id = i-1;	}
				}break;
				case 12:{
					patch_id = i;
				}break;
				case 13:{
					if (i%_Length <= _Length-2){ patch_id = i+1;	}
				}break;
				case 14:{
					if (i%_Length <= _Length-3){ patch_id = i+2;	}
				}break;
				case 15:{
					if ((i <= _PatchNum-_Length-1) && (i%_Length >= 2)){ patch_id = (i/_Length+1)*_Length + i%_Length-2;	}
				}break;
				case 16:{
					if ((i <= _PatchNum-_Length-1) && (i%_Length >= 1)){ patch_id = (i/_Length+1)*_Length + i%_Length-1;	}
				}break;
				case 17:{
					if (i <= _PatchNum-_Length-1){ patch_id = (i/_Length+1)*_Length + i%_Length;	}
				}break;
				case 18:{
					if ((i <= _PatchNum-_Length-1) && (i%_Length <= _Length-2)){ patch_id = (i/_Length+1)*_Length + i%_Length+1;	}
				}break;
				case 19:{
					if ((i <= _PatchNum-_Length-1) && (i%_Length <= _Length-3)){ patch_id = (i/_Length+1)*_Length + i%_Length+2;	}
				}break;
				case 20:{
					if ((i <= _PatchNum-_Length*2-1) && (i%_Length >= 2)){ patch_id = (i/_Length+2)*_Length + i%_Length-2;	}
				}break;
				case 21:{
					if ((i <= _PatchNum-_Length*2-1) && (i%_Length >= 1)){ patch_id = (i/_Length+2)*_Length + i%_Length-1;	}
				}break;
				case 22:{
					if (i <= _PatchNum-_Length*2-1){ patch_id = (i/_Length+2)*_Length + i%_Length;	}
				}break;
				case 23:{
					if ((i <= _PatchNum-_Length*2-1) && (i%_Length <= _Length-2)){ patch_id = (i/_Length+2)*_Length + i%_Length+1;	}
				}break;
				case 24:{
					if ((i <= _PatchNum-_Length*2-1) && (i%_Length <= _Length-3)){ patch_id = (i/_Length+2)*_Length + i%_Length+2;	}
				}break;
			}
			double offspring_trait = parent_trait;
			_Patch[patch_id].AddIndividual(patch_id, offspring_trait, 0);
		}
	}
}

void PopProcess::ErrorHandling1(){
	_TotalCoopNum   = 0;
	_OffspringNum	= 0;
	_CoopProp       = -1.0;
	_AverCoopDeg    = -1.0;	
}

void PopProcess::ErrorHandling2(){
    _PopSize        = -1;
	_TotalCoopNum   = -1;
	_OffspringNum	= -1;
	_CoopProp       = -1.0;
	_AverCoopDeg    = -1.0;
}

void PopProcess::ComputePopCoopData(){
	_TotalCoopNum = 0;
	_AverCoopDeg  = 0.0;
	for (Patch& patch: _Patch){
		patch.UpdateCoopNum();
		patch.UpdateCoopDegSum();
		_TotalCoopNum = _TotalCoopNum + patch.GetCoopNum();
		_AverCoopDeg  = _AverCoopDeg  + patch.GetCoopDegSum();
	}
	_CoopProp    = double(_TotalCoopNum)/double(_PopSize);
    _AverCoopDeg = _AverCoopDeg/double(_PopSize);
}

void PopProcess::ComputeTimeSeriesData(int Time){
	cout << _Replic << "," << Time+1 << "," << _HabitatAvail << "," << _Temperature << "," << _Length << "," << _ResourceAvail << "," << _InitPopDens << "," \
		 << _Topt << "," << _CTmax << "," << _CoopEfficiency << "," << _HalfConst << "," << _ReproMax << "," << _CostRate << "," << _Sociality << "," << _PopSize << "," \
		 << _OffspringNum << "," << _TotalCoopNum << "," << _PopSize-_TotalCoopNum << "," << _CoopProp << "," << _AverCoopDeg << "," << _InitCoopProp;
	
	for (int var_type = 0; var_type < 3; var_type++){
		for (int i = 0; i < _Patch.size(); i++){
			cout << ",";
			switch (var_type)
			{
				case 0:{
					cout << _PatchIsAvail[i];
				}break;
				case 1:{
					cout << _Patch[i].GetSize();
				}break;
				case 2:{
					cout << _Patch[i].GetOffspring().size();
				}break;
			}
		}
	}
	cout << "\n";
}

Simulation::Simulation(int OutputType, double HabitatAvail, double Temperature, int HabitatLength, double ResourceAvail, double InitPopDens, int EndTime, \
					   int InitReplic, int ReplicNum, int Span, double Topt, double CTmax, bool Sociality, double CoopEfficiency, double HalfConst, double ReproMax, \
					   double CostRate, double InitCoopProp, int RandomSeed){
	_OutputType				= OutputType;
	_HabitatAvail			= HabitatAvail;
	_Temperature			= Temperature;
	_HabitatLength			= HabitatLength;
	_ResourceAvail			= ResourceAvail;
	_InitPopDens			= InitPopDens;
	_EndTime				= EndTime;
	_InitReplic				= InitReplic;
	_ReplicNum				= ReplicNum;
	_Span					= Span;
	_Topt					= Topt;
	_CTmax					= CTmax;
	_Sociality			    = Sociality;
	_CoopEfficiency			= CoopEfficiency;
	_HalfConst				= HalfConst;
	_ReproMax				= ReproMax;
	_CostRate				= CostRate;
	_InitCoopProp			= InitCoopProp;
	_RandomSeed				= RandomSeed;
}

void Simulation::ProcessShunt(int OutputType){
	switch (OutputType)
	{
		case 1:{
			PopFinalState();
		}break;
		case 2:{
			PopTimeSeries();
		}break;
	}
}

void Simulation::PopFinalState(){
	for (int i = 0; i < _ReplicNum; i++){
		PopProcess PP(_OutputType, _HabitatAvail, _Temperature, _HabitatLength, _ResourceAvail, _InitPopDens, i, _Span, _Topt, _CTmax, _Sociality, _CoopEfficiency, \
					  _HalfConst, _ReproMax, _CostRate, _InitCoopProp, _RandomSeed);
		PP.SystemUpdate(_EndTime);
		
		cout << _InitReplic+i << "," << _HabitatAvail << "," << _Temperature << "," << _HabitatLength << "," << _ResourceAvail << "," << _InitPopDens << "," \
			 << _Topt << "," << _CTmax << "," << _CoopEfficiency << "," << _HalfConst << "," << _ReproMax << "," << _CostRate << "," << _Sociality << "," \
			 << PP.GetPopSize() << "," << PP.GetOffspringNum() << "," << PP.GetTotalCoopNum() << "," << PP.GetPopSize()-PP.GetTotalCoopNum() << "," \
			 << PP.GetCoopProp() << "," << PP.GetAverCoopDeg()  << "," << _InitCoopProp;

		vector<bool> patch_is_avail = PP.GetPatchIsAvail();
		vector<Patch> patches = PP.GetPatch();
		for (int var_type = 0; var_type < 3; var_type++){
			for (int i = 0; i < patches.size(); i++){
				cout << ",";
				switch (var_type)
				{
					case 0:{
						cout << patch_is_avail[i];
					}break;
					case 1:{
						cout << patches[i].GetSize();
					}break;
					case 2:{
						cout << patches[i].GetOffspring().size();
					}break;
				}
			}
		}
		cout << "\n";
	}
}

void Simulation::PopTimeSeries(){
	for (int i = 0; i < _ReplicNum; i++){
		PopProcess PP(_OutputType, _HabitatAvail, _Temperature, _HabitatLength, _ResourceAvail, _InitPopDens, i, _Span, _Topt, _CTmax, _Sociality, _CoopEfficiency, \
					  _HalfConst, _ReproMax, _CostRate, _InitCoopProp, _RandomSeed);
		
		cout << _InitReplic+i << ",0," << _HabitatAvail << "," << _Temperature << "," << _HabitatLength << "," << _ResourceAvail << "," << _InitPopDens << "," \
			 << _Topt << "," << _CTmax << "," << _CoopEfficiency << "," << _HalfConst << "," << _ReproMax << "," << _CostRate << "," << _Sociality << "," \
			 << PP.GetPopSize() << "," << PP.GetOffspringNum() << "," << PP.GetTotalCoopNum() << "," << PP.GetPopSize()-PP.GetTotalCoopNum() << "," \
			 << PP.GetCoopProp() << "," << PP.GetAverCoopDeg() << "," << _InitCoopProp;

		vector<bool> patch_is_avail = PP.GetPatchIsAvail();
		vector<Patch> patches = PP.GetPatch();
		for (int var_type = 0; var_type < 3; var_type++){
			for (int i = 0; i < patches.size(); i++){
				cout << ",";
				switch (var_type)
				{
					case 0:{
						cout << patch_is_avail[i];
					}break;
					case 1:{
						cout << patches[i].GetSize();
					}break;
					case 2:{
						cout << patches[i].GetOffspring().size();
					}break;
				}
			}
		}
		cout << "\n";

		PP.SystemUpdate(_EndTime);
	}
}