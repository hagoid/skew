#include <iostream>

using namespace std;

int rad(int n, int t){
	int a, i;
	for(i=1,a=1;(a = (a*t)%n)!=1;i++);
	return i;
}

int gcd(int n, int k){
	if(n%k == 0) return k;
	return gcd(k, n%k);
}

int lcm(int n, int k){
	return n*k / gcd(n,k);
}

class permutation{
	public:
		int n, cn;
		int *predpis;
		int *cyklus;
		bool *nastavene;
		int rad();
		void reset();
		void vypis();
		permutation(int tn);
};

permutation::permutation(int tn){
	n = tn;
	predpis = new int[n];
	nastavene = new bool[n];
	reset();
}

void permutation::reset(){
	for(int i=0; i<n; i++){
		predpis[i] = i;
		nastavene[i] = false;
	}
}

int permutation::rad(){
	int r = 1;
	int d, m;
	bool uz[n];
	for(int i=0; i<n; i++) uz[i] = false;
	for(int i=0; i<n; i++){
		if(!uz[i]){
			uz[i] = true;
			for(d=1,m=i;(m = predpis[m]) != i;d++) uz[m] = true;
			r = lcm(r, d);
		}
	}
	return r;
}

void permutation::vypis(){
	int m;
	bool uz[n], jeden;
	for(int i=0; i<n; i++) uz[i] = false;
	for(int i=0; i<n; i++){
		if(!uz[i]){
			m = i;
			jeden = true;
			while(predpis[m]!=i){
				if(jeden){
					cout<<"("<<i;
					uz[i] = true;
					jeden = false;
				}
				cout<<", ";
				m = predpis[m];
				uz[m] = true;
				cout << m;
			}
			if(!jeden) cout << ")";
		}
	}
	cout << "\n";
}

void permutuj(int k, int t, int h, permutation * phi){
	phi->reset();
	int i, p;
	for(i=0;i<phi->n;i+=k){
		phi->predpis[i] = (phi->predpis[i]*t) % phi->n;
		phi->nastavene[i] = true;
	}
	for(i=1, p=1; (phi->predpis[i] = (i + h)%phi->n) != 1; h=(h*t)%phi->n, i = phi->predpis[i], p++) phi->nastavene[i] = true;
	phi->nastavene[i] = true;
	phi->cn = p;
	phi->cyklus = new int[p];
	phi->cyklus[0] = 1;
	for(i=1;i<p;i++) phi->cyklus[i] = phi->predpis[phi->cyklus[i-1]];
}

bool over(int k, int t, int pi, permutation * phi){
	int mocnina = 1;
	int vysledok = phi->predpis[1];
	for(int i=1; i<k; i++){
		mocnina = (mocnina*pi)%phi->cn;
		vysledok = (vysledok + phi->cyklus[mocnina])%phi->n;
	}
	return (k*t)%phi->n == vysledok;
}

void dopermutuj(int pi, permutation * phi){
	int vysledok;
	int mocnina;
	for(int a=1; a<phi->n; a++) if(!phi->nastavene[a]){
		mocnina = 1;
		vysledok = phi->predpis[1];
		for(int i=1; i<a; i++){
			mocnina = (mocnina*pi)%phi->cn;
			vysledok = (vysledok + phi->cyklus[mocnina])%phi->n;
		}
		phi->predpis[a] = vysledok;
		phi->nastavene[a] = true;
	}
}

int main(){
	int n, k, t, h, r, pi;
	cin >> n;
	permutation phi(n);
	for(k = 2; k<=n; k++){
		if(n%k == 0){
			for(t = 1; t <= n/k; t++){
				if(gcd(n/k, t) == 1){
					for(h = 0; h < n; h += k){
						permutuj(k,t,h,&phi);
						r = phi.rad();
						for(pi = 1; pi<r; pi++){
							if(gcd(r, pi) == 1){
								if(rad(r, pi) == k){
									if(over(k,t,pi,&phi)){
										dopermutuj(pi,&phi);
										phi.vypis();
										cout <<n<<" "<<k<<" "<<t<<" "<<h<<" "<<pi<<"\n";
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
