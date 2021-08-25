
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "StringStreamUtils.h"

#include <cassert>
#include "RASCATGTRFiniteGammaPhyloProcess.h"
#include "Parallel.h"
#include <string>


void RASCATGTRFiniteGammaPhyloProcess::GlobalUpdateParameters()	{
	// MPI2
	// should send the slaves the relevant information
	// about model parameters

	// for this model, should broadcast
	// double alpha
	// int Ncomponent
	// int* alloc
	// double* rr
	// double** profile
	// double* brancharray
	// (but should first call PutBranchLengthsIntoArray())
	// 
	// upon receiving this information
	// slave should 
	// store it in the local copies of the variables
	// and then call
	// SetBranchLengthsFromArray()
	// SetAlpha(inalpha)

	assert(myid == 0);

	// ResampleWeights();
	RenormalizeProfiles();

	int i,j,nrr,nbranch = GetNbranch(),ni,nd,L1,L2;
	nrr = GetNrr();
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = 3 + nbranch + nrr + L2 + L1*(L2+1);
	ni = 1 + GetNsite();
	int ivector[ni];
	double dvector[nd]; 
	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	// GlobalBroadcastTree();

	// First we assemble the vector of doubles for distribution
	int index = 0;
	dvector[index] = GetAlpha();
	index ++;
	dvector[index] = branchalpha;
	index++;
	dvector[index] = branchbeta;
	index++;
	for(i=0; i<nbranch; ++i) {
		dvector[index] = blarray[i];
		index++;
	}
	
	for(i=0; i<nrr ; ++i) {
		dvector[index] = rr[i];
		index++;
	}

	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			dvector[index] = profile[i][j];
			index++;
		}
		dvector[index] = weight[i];
		index++;
	}
	for (int i=0; i<GetDim(); i++)	{
		dvector[index] = dirweight[i];
		index++;
	}

	// Now the vector of ints
	ivector[0] = GetNcomponent();
	for(i=0; i<GetNsite(); ++i) {
		ivector[1+i] = FiniteProfileProcess::alloc[i];
	}

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
}


void RASCATGTRFiniteGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	assert(myid > 0);

	switch(signal) {

	case UPDATE_RATE:
		SlaveUpdateRateSuffStat();
		break;
	case UPDATE_RRATE:
		SlaveUpdateRRSuffStat();
		break;
	case REALLOC_MOVE:
		SlaveIncrementalFiniteMove();
		break;
	case PROFILE_MOVE:
		SlaveMoveProfile();
		break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}

void RASCATGTRFiniteGammaPhyloProcess::SlaveUpdateParameters()	{

	// SlaveBroadcastTree();

	int i,j,L1,L2,ni,nd,nbranch = GetNbranch(),nrr = GetNrr();
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = 3 + nbranch + nrr + L2 + L1*(L2+1);
	ni = 1 + GetNsite();
	int* ivector = new int[ni];
	double* dvector = new double[nd];
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	int index = 0;
	SetAlpha(dvector[index]);
	index ++;
	branchalpha = dvector[index];
	index++;
	branchbeta = dvector[index];
	index++;
	for(i=0; i<nbranch; ++i) {
		blarray[i] = dvector[index];
		index++;
	}
	for(i=0; i<nrr; ++i) {
		rr[i] = dvector[index];
		index++;
	}
	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			profile[i][j] = dvector[index];
			index++;
		}
		weight[i] = dvector[index];
		index++;
	}
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = dvector[index];
		index++;
	}

	Ncomponent = ivector[0];
	for(i=0; i<GetNsite(); ++i) {
		FiniteProfileProcess::alloc[i] = ivector[1+i];
	}
	delete[] dvector;
	delete[] ivector;

	UpdateMatrices();
}

void RASCATGTRFiniteGammaPhyloProcess::ReadPB(int argc, char* argv[])	{

	string name = "";

	int burnin = -1;
	int every = 1;
	int until = -1;
	int ppred = 0;
	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic

	int cv = 0;
	int sitelogl = 0;
    int verbose = 0;
	int map = 0;
	int rates = 0;
	string testdatafile = "";
	int rateprior = 0;
	int profileprior = 0;
	int rootprior = 1;
	int savetrees = 0;

	int rr = 0;
	int ss = 0;

	int ancstatepostprobs = 0;

    int posthyper = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if (s == "-div")	{
				ppred = 2;
			}
			else if (s == "-comp")	{
				ppred = 3;
			}
            else if (s == "-siteconvprob")   {
                ppred = 4;
            }
            else if (s == "-sitecomp")   {
                ppred = 5;
            }
			else if (s == "-ppred")	{
				ppred = 1;
			}
			else if (s == "-allppred")	{
				ppred = -1;
			}
			else if (s == "-ppredrate")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					rateprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					rateprior = 0;
				}
				else	{
					cerr << "error after ppredrate: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-ppredprofile")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					profileprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					profileprior = 0;
				}
				else	{
					cerr << "error after ppredprofile: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-ppredroot")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					rootprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					rootprior = 0;
				}
				else	{
					cerr << "error after ppredroot: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-savetrees")	{
				savetrees = 1;
			}
			else if (s == "-anc")	{
				ancstatepostprobs = 1;
			}
			else if (s == "-ss")	{
				ss = 1;
			}
			else if (s == "-rr")	{
				rr = 1;
			}
			else if (s == "-sitelogl")	{
				sitelogl = 1;
			}
            else if (s == "-v") {
                verbose = 1;
            }
			else if (s == "-r")	{
				rates = 1;
			}
			else if (s == "-map")	{
				map = 1;
			}
            else if (s == "-posthyper") {
                posthyper = 1;
            }

			else if (s == "-oldcv")	{
				cv = 1;
				i++;
				testdatafile = argv[i];
			}

			else if ((s == "-cv") || (s == "-sitecv"))	{
				cv = 2;
				i++;
				testdatafile = argv[i];
			}

			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					string tmp = argv[i];
					if (IsInt(tmp))	{
						until = atoi(argv[i]);
					}
					else	{
						i--;
					}
				}
				else {
					i--;
				}
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
	}
	catch(...)	{
		cerr << "error in command\n";
		cerr << '\n';
		exit(1);
	}

	if (until == -1)	{
		until = GetSize();
	}
	if (burnin == -1)	{
		burnin = GetSize() / 5;
	}

	if ((GetNprocs() == 1) && (ppred || cv || sitelogl))	{
		cerr << "error : should run readpb_mpi in mpi mode, with at least 2 processes\n";
		MPI_Finalize();
		exit(1);
	}

	if (cv == 1)	{
		ReadCV(testdatafile,name,burnin,every,until);
	}
    else if (cv == 2)	{
		ReadSiteCV(testdatafile,name,burnin,every,until);
	}
	else if (sitelogl)	{
		ReadSiteLogL(name,burnin,every,until,verbose);
	}
	else if (ancstatepostprobs)	{
		ReadAncestral(name,burnin,every,until);
	}
	else if (ss)	{
		ReadSiteProfiles(name,burnin,every,until);
	}
	else if (rr)	{
		ReadRelRates(name,burnin,every,until);
	}
	else if (rates)	{
		ReadSiteRates(name,burnin,every,until);
	}
	else if (ppred == -1)	{
		AllPostPred(name,burnin,every,until,rateprior,profileprior,rootprior);
	}
	else if (ppred)	{
		PostPred(ppred,name,burnin,every,until,rateprior,profileprior,rootprior,savetrees);
	}
	else if (map)	{
		ReadMap(name,burnin,every,until);
	}
    else if (posthyper) {
		ReadPostHyper(name,burnin,every,until);
    }
	else	{
		Read(name,burnin,every,until);
	}
}

void RASCATGTRFiniteGammaPhyloProcess::SlaveComputeCVScore()	{

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}

	sitemax = sitemin + testsitemax - testsitemin;
	double** sitelogl = new double*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
		sitelogl[i] = new double[GetNcomponent()];
	}
	
	// UpdateMatrices();

	for (int k=0; k<GetNcomponent(); k++)	{
		for (int i=sitemin; i<sitemax; i++)	{
			ExpoConjugateGTRFiniteProfileProcess::alloc[i] = k;
		}
		UpdateConditionalLikelihoods();
		for (int i=sitemin; i<sitemax; i++)	{
			sitelogl[i][k] = sitelogL[i];
		}
	}

	double total = 0;
	for (int i=sitemin; i<sitemax; i++)	{
		double max = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			if ((!k) || (max < sitelogl[i][k]))	{
				max = sitelogl[i][k];
			}
		}
		double tot = 0;
		double totweight = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			tot += weight[k] * exp(sitelogl[i][k] - max);
			totweight += weight[k];
		}
		total += log(tot) + max;
	}

	MPI_Send(&total,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	for (int i=sitemin; i<sitemax; i++)	{
		delete[] sitelogl[i];
	}
	delete[] sitelogl;

	sitemax = bksitemax;

}

void RASCATGTRFiniteGammaPhyloProcess::SlaveComputeSiteLogL()	{

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}

	double** sitelogl = new double*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            sitelogl[i] = new double[GetNcomponent()];
        }
	}
	
	// UpdateMatrices();

	for (int k=0; k<GetNcomponent(); k++)	{
		for (int i=sitemin; i<sitemax; i++)	{
            if (ActiveSite(i))  {
                ExpoConjugateGTRFiniteProfileProcess::alloc[i] = k;
            }
		}
		UpdateConditionalLikelihoods();
		for (int i=sitemin; i<sitemax; i++)	{
            if (ActiveSite(i))  {
                sitelogl[i][k] = sitelogL[i];
            }
		}
	}

	double* meansitelogl = new double[GetNsite()];
    double* cumul = new double[GetNcomponent()];
	for (int i=0; i<GetNsite(); i++)	{
		meansitelogl[i] = 0;
	}
	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            double max = 0;
            for (int k=0; k<GetNcomponent(); k++)	{
                if ((!k) || (max < sitelogl[i][k]))	{
                    max = sitelogl[i][k];
                }
            }
            double tot = 0;
            for (int k=0; k<GetNcomponent(); k++)	{
                tot += weight[k] * exp(sitelogl[i][k] - max);
                cumul[k] = tot;
            }
            meansitelogl[i] = log(tot) + max;
            int k = 0;
            double u = tot*rnd::GetRandom().Uniform();
            while ((k < GetNcomponent()) && (u>cumul[k]))   {
                k++;
            }
            if (k == GetNcomponent())   {
                cerr << "error in RASCATFiniteGammaPhyloProcess::SlaveComputeSiteLogL: overflow\n";
                exit(1);
            }
            ExpoConjugateGTRFiniteProfileProcess::alloc[i] = k;
        }
	}
    UpdateConditionalLikelihoods();

	MPI_Send(meansitelogl,GetNsite(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	for (int i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            delete[] sitelogl[i];
        }
	}
	delete[] sitelogl;
	delete[] meansitelogl;
    delete[] cumul;
}


void RASCATGTRFiniteGammaPhyloProcess::ReadRelRates(string name, int burnin, int every, int until)	{

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

	double* meanrr = new double[GetNrr()];
	for (int k=0; k<GetNrr(); k++)	{
		meanrr[k] = 0;
	}

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		double total = 0;
		for (int k=0; k<GetNrr(); k++)	{
			total += rr[k];
		}
		total /= GetNrr();

		for (int k=0; k<GetNrr(); k++)	{
			meanrr[k] += rr[k] / total;
		}

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	for (int k=0; k<GetNrr(); k++)	{
		meanrr[k] /= samplesize;
	}
	ofstream os((name + ".meanrr").c_str());
	for (int k=0; k<GetDim(); k++)	{
		os << GetStateSpace()->GetState(k) << ' ';
	}
	os << '\n';
	os << '\n';
	for (int i=0; i<GetDim(); i++)	{
		for (int j=i+1; j<GetDim(); j++)	{
			os << GetStateSpace()->GetState(i) << '\t' << GetStateSpace()->GetState(j) << '\t' << meanrr[rrindex(i,j,GetDim())] << '\n';
		}
	}
	cerr << "mean relative exchangeabilities in " << name << ".meanrr\n";

}

void RASCATGTRFiniteGammaPhyloProcess::ReadSiteProfiles(string name, int burnin, int every, int until)	{

	double** sitestat = new double*[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		sitestat[i] = new double[GetDim()];
		for (int k=0; k<GetDim(); k++)	{
			sitestat[i][k] = 0;
		}
	}
	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		for (int i=0; i<GetNsite(); i++)	{
			double* p = GetProfile(i);
			for (int k=0; k<GetDim(); k++)	{
				sitestat[i][k] += p[k];
			}
		}
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	
	ofstream os((name + ".siteprofiles").c_str());
	for (int k=0; k<GetDim(); k++)	{
		os << GetStateSpace()->GetState(k) << ' ';
	}
	os << '\n';
	os << '\n';
	for (int i=0; i<GetNsite(); i++)	{
		os << i+1;
		for (int k=0; k<GetDim(); k++)	{
			sitestat[i][k] /= samplesize;
			os << '\t' << sitestat[i][k];
		}
		os << '\n';
	}
	cerr << "mean site-specific profiles in " << name << ".siteprofiles\n";
	cerr << '\n';
}

void RASCATGTRFiniteGammaPhyloProcess::ReadPostHyper(string name, int burnin, int every, int until)	{

    double meanratealpha = 0;
    double varratealpha = 0;
    vector<double> meandirweight(GetDim(), 0);
    vector<double> vardirweight(GetDim(), 0);
    vector<double> meanbl(GetNbranch(), 0);
    vector<double> varbl(GetNbranch(), 0);

    vector<double> meanfreq(GetDim(), 0);
    vector<double> varfreq(GetDim(), 0);

    vector<double> meanrr(GetNrr(), 0);
    vector<double> varrr(GetNrr(), 0);

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

        meanratealpha += alpha;
        varratealpha += alpha*alpha;
        for (int k=0; k<GetDim(); k++)  {
            meandirweight[k] += dirweight[k];
            vardirweight[k] += dirweight[k]*dirweight[k];
        }

        if (fixncomp && (GetNcomponent() == 1)) {
            for (int k=0; k<GetDim(); k++)  {
                meanfreq[k] += profile[0][k];
                varfreq[k] += profile[0][k] * profile[0][k];
            }
        }

        // double meanRR = 1.0;
        if (! fixrr)    {
            /*
            for (int k=0; k<GetNrr(); k++)  {
                meanRR += rr[k];
            }
            meanRR /= GetNrr();
            */
            for (int k=0; k<GetNrr(); k++)  {
                double tmp = rr[k];
                // double tmp = rr[k] / meanRR;
                meanrr[k] += tmp;
                varrr[k] += tmp*tmp;
            }
        }

        for (int j=1; j<GetNbranch(); j++)  {
            double tmp = blarray[j];
            // double tmp = blarray[j] * meanRR;
            meanbl[j] += tmp;
            varbl[j] += tmp*tmp;
            // meanbl[j] += blarray[j];
            // varbl[j] += blarray[j]*blarray[j];
        }

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	
	ofstream os((name + ".posthyper").c_str());
    meanratealpha /= samplesize;
    varratealpha /= samplesize;
    varratealpha -= meanratealpha*meanratealpha;
    os << meanratealpha*meanratealpha / varratealpha << '\t';
    os << meanratealpha / varratealpha << '\n';
    os << '\n';
	if (! dirweightprior)   {
        for (int k=0; k<GetDim(); k++)  {
            meandirweight[k] /= samplesize;
            vardirweight[k] /= samplesize;
            vardirweight[k] -= meandirweight[k]*meandirweight[k];
            os << meandirweight[k]*meandirweight[k]/vardirweight[k] << '\t';
            os << meandirweight[k]/vardirweight[k] << '\n';
        }
        os << '\n';
    }
    if (fixncomp && (GetNcomponent() == 1)) {
        double num = 0;
        double var = 0;
        for (int k=0; k<GetDim(); k++)  {
            meanfreq[k] /= samplesize;
            num += meanfreq[k] * (1-meanfreq[k]);
            varfreq[k] /= samplesize;
            varfreq[k] -= meanfreq[k]*meanfreq[k];
            var += varfreq[k];
        }
        double conc = num/var - 1;
        for (int k=0; k<GetDim(); k++)  {
            os << conc * meanfreq[k] << '\t';
        }
        os << '\n';
        os << '\n';
    }

    if (! fixrr)    {
        for (int k=0; k<GetNrr(); k++)  {
            meanrr[k] /= samplesize;
            varrr[k] /= samplesize;
            varrr[k] -= meanrr[k]*meanrr[k];
            os << meanrr[k]*meanrr[k]/varrr[k] << '\t';
            os << meanrr[k]/varrr[k] << '\n';
        }
        os << '\n';
    }

    for (int j=1; j<GetNbranch(); j++)  {
        meanbl[j] /= samplesize;
        varbl[j] /= samplesize;
        varbl[j] -= meanbl[j]*meanbl[j];
        os << meanbl[j]*meanbl[j]/varbl[j] << '\t';
        os << meanbl[j]/varbl[j] << '\n';
    }
	cerr << "posterior mean shape and scale params in " << name << ".posthyper\n";
	cerr << '\n';
}

void RASCATGTRFiniteGammaPhyloProcess::GlobalSetEmpiricalPrior(istream& is)    {

    // read from stream
    is >> empalpha >> empbeta;
    if (! dirweightprior)   {
        for (int k=0; k<GetDim(); k++)  {
            is >> empdirweightalpha[k] >> empdirweightbeta[k];
        }
    }
    if (fixncomp && (GetNcomponent() == 1)) {
        for (int k=0; k<GetDim(); k++)  {
            is >> empdirweight[k];
        }
    }
    if (! fixrr)    {
        for (int k=0; k<GetNrr(); k++)  {
            is >> emprralpha[k] >> emprrbeta[k];
        }
        /*
        double pseudocount = 1.0;
        for (int k=0; k<GetNrr(); k++)  {
            emprralpha[k] += pseudocount * LG_RR[k];
            emprrbeta[k] += pseudocount * LG_RR[k];
        }
        */
    }
    for (int j=1; j<GetNbranch(); j++)  {
        is >> branchempalpha[j] >> branchempbeta[j];
    }
	MESSAGE signal = EMPIRICALPRIOR;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&empalpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&empbeta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (! dirweightprior)   {
        MPI_Bcast(empdirweightalpha,GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(empdirweightbeta,GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    if (fixncomp && (GetNcomponent() == 1)) {
        MPI_Bcast(empdirweight,GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    if (! fixrr)    {
        MPI_Bcast(emprralpha,GetNrr(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(emprrbeta,GetNrr(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
	MPI_Bcast(branchempalpha,GetNbranch(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(branchempbeta,GetNbranch(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void RASCATGTRFiniteGammaPhyloProcess::SlaveSetEmpiricalPrior()    {

	MPI_Bcast(&empalpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&empbeta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (! dirweightprior)   {
        MPI_Bcast(empdirweightalpha,GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(empdirweightbeta,GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    if (fixncomp && (GetNcomponent() == 1)) {
        MPI_Bcast(empdirweight,GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    if (! fixrr)    {
        MPI_Bcast(emprralpha,GetNrr(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(emprrbeta,GetNrr(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
	MPI_Bcast(branchempalpha,GetNbranch(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(branchempbeta,GetNbranch(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

double RASCATGTRFiniteGammaPhyloProcess::GlobalGetSiteSteppingLogLikelihoodNonIS(int site, int nrep0, int restore) {

    int nrep_per_proc = nrep0 / (GetNprocs()-1);
    if (nrep0 % (GetNprocs() - 1))  {
        nrep_per_proc ++;
    }
    // int nrep = nrep_per_proc * (GetNprocs()-1);

    MESSAGE signal = STEPPINGSITELOGL;
    MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
    int param[3];
    param[0] = site;
    param[1] = nrep_per_proc;
    param[2] = restore;
    MPI_Bcast(param,3,MPI_INT,0,MPI_COMM_WORLD);

    int oldalloc = ExpoConjugateGTRFiniteProfileProcess::alloc[site];

    double master_logl[GetNprocs()];
    double slave_logl;

    int master_alloc[GetNprocs()];
    int slave_alloc = -1;

    MPI_Gather(&slave_logl, 1, MPI_DOUBLE, master_logl, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&slave_alloc, 1, MPI_INT, master_alloc, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double max = 0;
    for (int i=1; i<GetNprocs(); i++)   {
        if (master_alloc[i] != -1)    {
            if ((!max) || (max < master_logl[i]))    {
                max = master_logl[i];
            }
        }
    }
    double tot = 0;
    double post[GetNprocs()-1];
    for (int i=1; i<GetNprocs(); i++)   {
        if (master_alloc[i] != -1)    {
            double tmp = exp(master_logl[i]-max);
            post[i-1] = tmp;
            tot += tmp;
        }
        else    {
            post[i-1] = 0;
        }
    }
    if (! tot)  {
        cerr << "error in stepping logl: total likelihood is 0\n";
        exit(1);
    }

    double L = log(tot) + max;

    for (int i=1; i<GetNprocs(); i++)   {
        post[i-1] /= tot;
    }
    int procalloc = rnd::GetRandom().FiniteDiscrete(GetNprocs()-1, post) + 1;

    if (! restore)  {
        int newalloc = master_alloc[procalloc];
        RemoveSite(site, oldalloc);
        AddSite(site, newalloc);
        GlobalUpdateParameters();
    }
    return L;
}

void RASCATGTRFiniteGammaPhyloProcess::SlaveGetSiteSteppingLogLikelihoodNonIS()    {

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}

    int param[3];
    MPI_Bcast(param,3,MPI_INT,0,MPI_COMM_WORLD);
    int site = param[0];
    // int nrep = param[1];
    int restore = param[2];

    int bkalloc = ExpoConjugateGTRFiniteProfileProcess::alloc[site];

	int width = GetNcomponent() / (GetNprocs()-1);
    int r = GetNcomponent() % (GetNprocs()-1);
	int smin[GetNprocs()-1];
	int smax[GetNprocs()-1];
    int s = 0;
	for(int i=0; i<GetNprocs()-1; i++) {
		smin[i] = s;
        if (i < r)  {
            s += width + 1;
        }
        else    {
            s += width;
        }
        smax[i] = s;
    }
    if (s != GetNcomponent())  {
        cerr << "error: ncomp checksum\n";
        exit(1);
    }

    int kmin = smin[myid-1];
    int kmax = smax[myid-1];
    int krange = kmax - kmin;

    double sitelogl[krange];
    for (int k=kmin; k<kmax; k++)	{
        ExpoConjugateGTRFiniteProfileProcess::alloc[site] = k;
        double tmp = SiteLogLikelihood(site);
        sitelogl[k-kmin] = tmp;
    }

    double max = 0;
    for (int k=kmin; k<kmax; k++)	{
        if ((!max) || (max < sitelogl[k-kmin]))	{
            max = sitelogl[k-kmin];
        }
    }

    double post[krange];
    double tot= 0;
    for (int k=kmin; k<kmax; k++)   {
        post[k-kmin] = weight[k] * exp(sitelogl[k-kmin] - max);
        tot += post[k-kmin];
    }

    double slave_logl = 0;
    if (tot > 0)    {
        slave_logl = log(tot) + max;
    }

    int slave_alloc = -1;
    if (tot > 0)    {
        for (int k=kmin; k<kmax; k++)   {
            post[k-kmin] /= tot;
        }
        slave_alloc = rnd::GetRandom().FiniteDiscrete(krange, post) + kmin;
    }

    double master_logl[GetNprocs()];
    int master_alloc[GetNprocs()];

    MPI_Gather(&slave_logl, 1, MPI_DOUBLE, master_logl, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&slave_alloc, 1, MPI_INT, master_alloc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (restore)    {
        ExpoConjugateGTRFiniteProfileProcess::alloc[site] = bkalloc;
    }
}

