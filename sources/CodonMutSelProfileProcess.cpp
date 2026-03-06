
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "CodonMutSelProfileProcess.h"
// include things here...


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//	* CodonMutSelProfileProcess
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

void CodonMutSelProfileProcess::Create(int innsite, int indim, CodonStateSpace* instatespace)	{
	
	if ( (!nucrr) && (!nucstat) && (!omega) )	{
		ProfileProcess::Create(innsite, indim);
		Nnucrr = Nnuc * (Nnuc-1) / 2;
		nucrr = new double[Nnucrr];
		nucstat = new double[Nnuc];
		statespace = instatespace;
		omega = new double;
		SampleNucRR();
		SampleNucStat();
		SampleOmega();				

	}
	else	{
		cerr << "Create of CodonMutSelProfileProcess, nucrr and/or nucstat are/is not 0.\n";
		exit(1);
	}
}


void CodonMutSelProfileProcess::Delete()	{

	if ( (nucrr) && (nucstat) && (omega) )	{
		delete[] nucrr;
		delete[] nucstat;
		delete omega;
		nucrr = 0;
		nucstat = 0;
		omega = 0;
		ProfileProcess::Delete();
	}
	else	{
		cerr << "Delete of CodonMutSelProfileProcess, nucrr and/or nucstat and/or omega are/is 0.\n";
		exit(1);
	}
}


double CodonMutSelProfileProcess::LogNucRRPrior()	{
	double total = 0;
	// exp prior
	//for (int i=0; i<GetNnucrr(); i++)  {
	//	total -= nucrr[i];
	//}
	return total;
}


void CodonMutSelProfileProcess::SampleNucRR()	{
	//cerr << "nucrr set for testing...\n";
	double total = 0;
	for (int i=0; i<GetNnucrr(); i++)  {
		//nucrr[i] = rnd::GetRandom().sExpo();
		// nucrr[i] = Random::sExpo();
		//*
		if (i==0) nucrr[i] = 1.0; // AC
		if (i==1) nucrr[i] = 1.0; // AG
		if (i==2) nucrr[i] = 1.0; // AT
		if (i==3) nucrr[i] = 1.0; // CG
		if (i==4) nucrr[i] = 1.0; // CT
		if (i==5) nucrr[i] = 1.0; // GT
		//*/
		total += nucrr[i];
	}
	for (int i=0; i<GetNnucrr(); i++)  {
		nucrr[i] /= total;
	}
}


double CodonMutSelProfileProcess::LogNucStatPrior()	{
	return 0;
}


void CodonMutSelProfileProcess::SampleNucStat()	{

	//cerr << "nucstat set to specific values for testing...\n";
	double total = 0;
	for (int i=0; i<Nnuc; i++)  {
		//nucstat[i] = rnd::GetRandom().sExpo();
		// nucstat[i] = Random::sExpo();
		//*
		if (i==0) nucstat[i] = 0.25;
		if (i==1) nucstat[i] = 0.25;
		if (i==2) nucstat[i] = 0.25;
		if (i==3) nucstat[i] = 0.25;
		//*/
		total += nucstat[i];
	}
	for (int i=0; i<Nnuc; i++)  {
		nucstat[i] /= total;
	}
}

double CodonMutSelProfileProcess::LogOmegaPrior()        {

	// flat
	//return 0;

	// exponential
	// return - *omega;
	
	// jeffreys
	// return -log(*omega);

	if (omegaprior == 0)	{
		// ratio of exponential random variables
		return -2 * log(1 + *omega);
	}
    else if (omegaprior == 1)   {
		// jeffreys
		return -log(*omega);
	}
    else if (omegaprior == 2)   {
        // gamma
        return -*omega;
    }
    else    {
        cerr << "error: did not recognize omega prior\n";
        exit(1);
    }
}

void CodonMutSelProfileProcess::SampleOmega()        {

	//double rv1 = -log(Random::Uniform());
	//double rv2 = -log(Random::Uniform());
	//*omega = rv1/rv2;
	*omega = 1.0;
}

double CodonMutSelProfileProcess::MoveOmega(double tuning)        {

	int naccepted = 0;
	double bkomega = *omega;
	double deltalogprob = -ProfileSuffStatLogProb();
	//double deltalogprob = 0;
	deltalogprob -= LogOmegaPrior();

	double h = tuning * (rnd::GetRandom().Uniform() -0.5);
	double e = exp(h);
	*omega *= e;


	UpdateMatrices();
	deltalogprob += h;
	deltalogprob += LogOmegaPrior();
	//cerr << "before calling ProfileSuffStatLogProb()\t";
	//cerr.flush();

	deltalogprob += ProfileSuffStatLogProb();

	//cerr << "in move, omega is " << *omega << "\n";
	//cerr.flush();
	//cerr << "deltalobprob is : " << deltalogprob << "\n";
	//cerr.flush();

	int accepted = (rnd::GetRandom().Uniform() < exp(deltalogprob));
	if (accepted)	{
		naccepted++;	
	}
	else	{
		*omega = bkomega;
		UpdateMatrices();
	}
	return naccepted;	
}

double CodonMutSelProfileProcess::MoveNucRR(double tuning)	{

	//exponental prior

	int naccepted = 0;
        for (int i=0; i<GetNnucrr(); i++)  {
                double deltalogratio = - LogNucRRPrior() - ProfileSuffStatLogProb();
                double bk = nucrr[i];
                double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
                // double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                nucrr[i] *= e;
                UpdateMatrices();
                deltalogratio += LogNucRRPrior() + ProfileSuffStatLogProb();
                deltalogratio += m;
                int accepted = (rnd::GetRandom().Uniform() < exp(deltalogratio));
                // int accepted = (Random::Uniform() < exp(deltalogratio));
                if (accepted)   {
                        naccepted++;
                }
                else    {
                        nucrr[i] = bk;
                        UpdateMatrices();
                }
        }
        return ((double) naccepted) / GetNnucrr();
}

double CodonMutSelProfileProcess::MoveNucRR(double tuning, int n)	{

	//dirichlet prior

	int naccepted = 0;
	double* bk = new double[GetNnucrr()];
	for (int k=0; k<GetNnucrr(); k++)  {
		bk[k] = nucrr[k];
	}
	double deltalogprob = -ProfileSuffStatLogProb();
	double loghastings = ProfileProposeMove(nucrr,tuning,n,GetNnucrr());
	UpdateMatrices();
	deltalogprob += ProfileSuffStatLogProb();
	deltalogprob += loghastings;
	int accepted = (rnd::GetRandom().Uniform() < exp(deltalogprob));
	// int accepted = (Random::Uniform() < exp(deltalogprob));
	if (accepted)   {
		naccepted ++;
	}
	else    {
		for (int k=0; k<GetNnucrr(); k++)  {
			nucrr[k] = bk[k];
		}
		UpdateMatrices();
	}
	delete[] bk;
	return naccepted; 
}


double CodonMutSelProfileProcess::MoveNucStat(double tuning, int n)	{

	int naccepted = 0;
	double* bk = new double[Nnuc];
	for (int k=0; k<Nnuc; k++)  {
		bk[k] = nucstat[k];
	}
	double deltalogprob = -ProfileSuffStatLogProb();
	double loghastings = ProfileProposeMove(nucstat,tuning,n,Nnuc);
	UpdateMatrices();
	deltalogprob += ProfileSuffStatLogProb();
	deltalogprob += loghastings;
	int accepted = (rnd::GetRandom().Uniform() < exp(deltalogprob));
	// int accepted = (Random::Uniform() < exp(deltalogprob));
	if (accepted)   {
		naccepted ++;
	}
	else    {
		for (int k=0; k<Nnuc; k++)  {
			nucstat[k] = bk[k];
		}
		UpdateMatrices();
	}
	delete[] bk;
	return naccepted; 
}


