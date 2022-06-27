#include <iostream>
#include <fstream>
#include <cstdio>

using namespace std;

class orbit{
	public:
		int * orbita;
		int * pi;
		int dlzka;
		orbit(int tdlzka);
		~orbit();
};

class permutation{
	public:
		int n;
		orbit ** orbity;
		int * poradie;
		int nasledovnik(int a, int k);
		void reset();
		void nacitaj();
		void vypis();
		void vypis_pi();
		void pi_na_orbite();
		void dopocitaj_z_orbity();
		void z_predpisu(int * predpis, int * pi);
		permutation(int tn);
		~permutation();
};

orbit::orbit(int tdlzka){
	dlzka = tdlzka;
	orbita = new int[dlzka];
	pi = new int[dlzka];
	for(int i=0; i<dlzka; i++){
		pi[i] = 0;
	}
}

orbit::~orbit(){
	delete orbita;
	delete pi;
}

permutation::permutation(int tn){
	n = tn;
	orbity = new orbit *[n];
	poradie = new int[n];
	for(int i=0; i<n; i++){
		orbity[i] = NULL;
		poradie[i] = 0;
	}
}

permutation::~permutation(){
	reset();
	delete poradie;
	delete orbity;
}

void permutation::reset(){
	bool znicene[n];
	int i;
	for(i=0; i<n; i++){
		if(poradie[i] == 0) delete orbity[i];
		orbity[i] = NULL;
		poradie[i] = 0;
	}
}

void permutation::z_predpisu(int * predpis, int * pi){
	int dlzka;
	for(i=0; i<n; i++){
		if(orbity[i] == NULL){
			dlzka = 0;
			for() OTESTUJ TO
		}
	}
}

int permutation::nasledovnik(int a, int k = 1){
	if(orbity[a] == NULL) return a;
	return orbity[a]->orbita[(poradie[a]+k)%orbity[a]->dlzka];
}

void permutation::nacitaj(){
	int p, d;
	orbit * o;
	cin >> p;
	for(int i=0; i<p; i++){
		cin >> d;
		o = new orbit(d);
		for(int j=0; j<d; j++){
			cin >> o->orbita[j];
			orbity[ o->orbita[j] ] = o;
			poradie[ o->orbita[j] ] = j;
		}
	}
}

void permutation::vypis(){
	bool nevypisane[n];
	for(int i=0; i<n; i++){
		if(orbity[i] && !poradie[i]){
			cout << "(" << orbity[i]->orbita[0];
			for(int j=1; j<orbity[i]->dlzka; j++){
				cout << ", " << orbity[i]->orbita[j];
			}
			cout << ")";
		}
	}
}

void permutation::vypis_pi(){
	bool nevypisane[n];
	for(int i=0; i<n; i++){
		if(orbity[i] && !poradie[i]){
			cout << "[" << orbity[i]->pi[0];
			for(int j=1; j<orbity[i]->dlzka; j++){
				cout << ", " << orbity[i]->pi[j];
			}
			cout << "]";
		}
	}
}

void permutation::pi_na_orbite(){
	if(orbity[1] == NULL || orbity[n-nasledovnik(1)] == NULL){
		printf("Neni orbita jednotky.\n");
		return;
	}
	int j;
	for(int i=0; i < orbity[1]->dlzka; i++){
		if(orbity[n-nasledovnik(1)] != orbity[n-1]){
			printf("Neni su orbity spravne.\n");
			return ;
		}
		j = orbity[1]->orbita[i];
		orbity[1]->pi[i] = (poradie[n-nasledovnik(j)] - poradie[n-j] + orbity[n-j]->dlzka) % orbity[n-j]->dlzka;
	}
}

void permutation::dopocitaj_z_orbity(){
	if(orbity[1] == NULL){
		//printf("Neni orbita jednotky.\n");
		return;
	}
	int d = orbity[1]->dlzka;
	int pi[n], predpis[n], sucet[d+1], i;
	sucet[0] = 0;
	for(i=0; i<d; i++){
		sucet[i+1] = sucet[i] + orbity[1]->pi[(poradie[1]+i)%d];
	}
	pi[1] = sucet[1];
	for(i=2; i<n; i++){
		pi[i] = (sucet[pi[i-1]%d] + sucet[d]*(pi[i-1]/d)) % d;
	}
	pi[0] = (sucet[pi[n-1]%d] + sucet[d]*(pi[n-1]/d)) % d;
	if(pi[0] != 1){
		//printf("Nesedi nam pi(0).\n");
		return;
	}
	for(i=1; i<d; i++){
		if(pi[orbity[1]->orbita[i]] != orbity[1]->pi[i]){
			//printf("Nesedia nam pi.\n");
			return;
		}
	}
	
	predpis[1] = nasledovnik(1);
	for(i=2; i<n; i++){
		predpis[i] = (predpis[i-1] + nasledovnik(1, pi[i-1])) % n;
	}
	predpis[0] = (predpis[n-1] + nasledovnik(1, pi[n-1])) % n;
	if(predpis[0] != 0){
		//printf("Nesedi nam predpis 0.\n");
		return;
	}
	for(i=0; i<n; i++){
		if(orbity[i] != NULL){
			if(predpis[i] != nasledovnik(i)){
				//printf("Nesedia nam predpisy.\n");
				return;
			}
		}
	}
	
	reset();
	z_predpisu(predpis, pi);
	
	//printf("Vsetko ok.\n");
}

int main(){
	int n, pocet;
	cin >> n >> pocet;
	permutation * phi = new permutation(n);
	for(int i=0; i<pocet; i++){
		phi->reset();
		phi->nacitaj();
		phi->pi_na_orbite();
		phi->dopocitaj_z_orbity();
	}
	return 0;
}

