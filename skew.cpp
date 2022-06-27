#include <iostream>
#include <fstream>
#include <cstdio>

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
		int *pi;
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
	pi = new int[n];
	cyklus = NULL;
	reset();
}

void permutation::reset(){
	for(int i=0; i<n; i++){
		predpis[i] = i;
		nastavene[i] = false;
		pi[i] = 1;
	}
	delete cyklus;
	cyklus = NULL;
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

int permutuj(int k, int t, int h, permutation * phi){
	phi->reset();
	int i, p, rH;
	for(i=0;i<phi->n;i+=k){
		phi->predpis[i] = (phi->predpis[i]*t) % phi->n;
		phi->nastavene[i] = true;
	}
	rH = phi->rad();
	for(i=1, p=1; (phi->predpis[i] = (i + h)%phi->n) != 1; h=(h*t)%phi->n, i = phi->predpis[i], p++) phi->nastavene[i] = true;
	phi->nastavene[i] = true;
	phi->cn = p;
	phi->cyklus = new int[p];
	phi->cyklus[0] = 1;
	for(i=1;i<p;i++) phi->cyklus[i] = phi->predpis[phi->cyklus[i-1]];
	return rH;
}

bool over(int k, int t, int pi, int rH, permutation * phi){
	int mocnina = 1;
	int vysledok = phi->predpis[1];
	for(int i=1; i<k; i++){
		mocnina = (mocnina*pi)%phi->cn;
		vysledok = (vysledok + phi->cyklus[mocnina])%phi->n;
	}
	return ((k*t)%phi->n == vysledok) && ((pi-1)%rH == 0);
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
	phi->pi[1] = pi;
}

int delitelov(int n){
	int pocet = 0;
	for(int i=1;i<=n;i++) if(n%i==0) pocet++;
	return pocet;
}

bool porovnaj(permutation *a, permutation *b){
	bool ok = true;
	if(a->n != b->n) return false;
	for(int i = 0;i<a->n;i++){
		ok &= (a->predpis[i] == b->predpis[i]);
	}
	ok &= (a->pi[1] == b->pi[1]);
	return ok;
}

int najdi(int n, int pocet = 0, permutation **zoznam = NULL){
	int k, t, h, r, rH, pi, vysledok=0;
	bool check[pocet], checked;
	for(int i=0;i<pocet;i++) check[i]=false;
	permutation *phi = new permutation(n);
	//cout << "Nasli sme:\n";
	for(k = 2; k<=n; k++){
		if(n%k == 0){
			for(t = 1; t <= n/k; t++){
				if(gcd(n/k, t) == 1){
					for(h = 0; h < n; h += k){
						rH = permutuj(k,t,h,phi);
						r = phi->rad();
						for(pi = 1; pi<r; pi++){
							if(gcd(r, pi) == 1){
								if(rad(r, pi) == k){
									if(over(k,t,pi,rH,phi)){
										dopermutuj(pi,phi);
										if(pocet){
											checked = false;
											for(int i=0;i<pocet;i++){
												if(porovnaj(zoznam[i],phi)){
													check[i] = true;
													checked = true;
												}
											}
											if(!checked){
												phi->vypis();
												cout <<n<<" "<<k<<" "<<t<<" "<<h<<" "<<pi<<"\n";
											}
										}
										else if(phi->rad()!=n/k){
											
											//phi->vypis();
											//cout <<n<<" "<<k<<" "<<t<<" "<<h<<" "<<pi<<"\n";
										}
										vysledok++;
										phi->reset();
										permutuj(k,t,h,phi);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//cout << "Nenasli sme:\n";
	for(int i=0;i<pocet;i++) if(!check[i]) zoznam[i]->vypis();
	return vysledok;
}

bool coel(permutation *phi){
	bool ok = true;
	ok &= (phi->pi[1]!=1);
	for(int i=1; i<phi->cn; i++) ok &= (phi->pi[1] == (phi->pi[phi->cyklus[i]]));
	return ok;
}

int main(){
	ifstream f;
	f.open("zoznam/zoznam");
	int n, pocet, cyklov, dlzka, naozaj, a, b , z;
	permutation *zoznam[100];
	permutation *phi;
	/*for(int i = 2; i<=52; i++){
		naozaj = 0;
		f >> n >> pocet;
		cout << n << "\n";
		phi = new permutation(n);
		for(int j=0; j<pocet; j++){
			phi->reset();
			f >> cyklov;
			if(cyklov == 0){
				phi->cn = 1;
				phi->cyklus = new int[1];
				phi->cyklus[0] = 1;
			}
			for(int k=0; k<cyklov; k++){
				f >> dlzka;
				f >> z;
				if(k == 0){
					phi->cn = dlzka;
					phi->cyklus = new int[dlzka];
					phi->cyklus[0] = 1;
				}
				a = z;
				for(int l=1; l<dlzka; l++){
					f >> b;
					if(k == 0) phi->cyklus[l] = b;
					phi->predpis[a] = b;
					a = b;
				}
				phi->predpis[a] = z;
			}
			for(int k=0; k<n; k++){
				f >> phi->pi[k];
			}
			if(coel(phi)){
				zoznam[naozaj++] = phi;
				phi = new permutation(n);
			}
		}
		
		najdi(n, naozaj, zoznam);
		for(int j=0; j<naozaj; j++) delete zoznam[j];
	}*/
	f.close();
	
	for(n=2;n<=300;n++){
		pocet = najdi(n);
		if(pocet) cout << "Pre n="<<n<<" je pocet " << pocet<<"\n";
	}
}
