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
		int najdi_k();
		void reset();
		void vypis(bool, int);
		void vypispower(bool, int);
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

int permutation::najdi_k(){
	int r = 1;
	int d, m, p;
	bool uz[n], opak;
	for(int i=0; i<n; i++) uz[i] = false;
	for(int i=0; i<n; i++){
		if(!uz[i]){
			uz[i] = true;
			m = i;
			d = 1;
			p = pi[i];
			opak = false;
			while((m = predpis[m]) != i){
				if(!opak){
					if(pi[m] == p) opak = true;
					else d++;
				}
				uz[m] = true;
			}
			r = lcm(r, d);
		}
	}
	return r;
}

void permutation::vypis(bool enter = true, int k = 1){
	int m;
	bool uz[n], jeden;
	for(int i=0; i<n; i++) uz[i] = false;
	for(int i=0; i<n; i++){
		if(!uz[i]){
			m = i;
			for(int j = 0; j < k-1; j++) m = predpis[m];
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
				for(int j = 0; j < k-1; j++) m = predpis[m];
			}
			if(!jeden) cout << ")";
		}
	}
	if(enter) cout << "\n";
}

void permutation::vypispower(bool enter = true, int k = 1){
	int m, p;
	bool uz[n], jeden;
	for(int i=0; i<n; i++) uz[i] = false;
	for(int i=0; i<n; i++){
		if(!uz[i]){
			m = i;
			p = pi[m];
			for(int j = 1; j < k; j++){
				m = predpis[m];
				p = p + pi[m];
			}
			jeden = true;
			while(predpis[m]!=i){
				if(jeden){
					cout<<"["<< (p/k)%(rad()/k);
					uz[i] = true;
					jeden = false;
				}
				cout<<", ";
				m = predpis[m];
				uz[m] = true;
				p = pi[m];
				for(int j = 1; j < k; j++){
					m = predpis[m];
					p = p + pi[m];
				}
				cout << (p/k)%(rad()/k);
			}
			if(!jeden) cout << "]";
		}
	}
	if(enter) cout << "\n";
}

int permutuj(int k, int t, int h, permutation * phi){
	phi->reset();
	int i, p, rH;
	for(i=0;i<phi->n;i+=k){
		phi->predpis[i] = (phi->predpis[i]*t) % phi->n;
		phi->nastavene[i] = true;
	}
	rH = phi->rad();
	if(k != 1){
		for(i=1, p=1; (phi->predpis[i] = (i + h)%phi->n) != 1; h=(h*t)%phi->n, i = phi->predpis[i], p++) phi->nastavene[i] = true;
		phi->nastavene[i] = true;
	}
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
	int r = phi->rad();
	for(int i = 2; i < phi->n; i++) phi->pi[i] = (phi->pi[i-1] * phi->pi[1])%r;
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

int * orbituj(permutation *phi){
	int * vysledok = new int[phi->n];
	int i, t, p;
	for(i = 0; i < phi->n; i++) vysledok[i] = 0;
	for(i = 0; i < phi->n; i++){
		if(vysledok[i] == 0){
			p = 1;
			t = phi->predpis[i];
			while(t != i){
				p++;
				t = phi->predpis[t];
			}
			t = phi->predpis[i];
			vysledok[i] = p;
			while(t != i){
				vysledok[t] = p;
				t = phi->predpis[t];
			}
			
		}
	}
	return vysledok;
}

void cekujorbity(int k, permutation * phi){
	int n = phi->n;
	int * orbity = orbituj(phi);
	int ref;
	for(int i = 2; i < k; i++){
		ref = orbity[i];
		for(int j = 1; j < n/k; j++){
			if(orbity[j*k+i] != ref){
				cout << n << " " << k <<"\n";
				phi->vypis();
			}
		}
	}
	delete orbity;
}

void peknevypis(int k, int t, permutation * phi){
	printf("C.O.P.F. ");
	if(k==1 && t==1) printf("Id(sn)\n");
	else phi->vypis();
	printf("  of Z_%d,\n", phi->n);
	printf("  with power function values ");
	phi->vypispower();
	printf("  and parameters k=%d, t=%d\n\n",k,t);
}

int najdi(int n, int pocet = 0, permutation **zoznam = NULL){
	int k, t, h, r, rH, pi, vysledok=0, ncopf=0, nauto=0;
	bool check[pocet], checked;
	for(int i=0;i<pocet;i++) check[i]=false;
	permutation *phi = new permutation(n);
	//cout << "Nasli sme:\n";
	for(k = 1; k<=n; k++){
		if(n%k == 0){
			for(t = 1; t <= n/k; t++){
				if(gcd(n/k, t) == 1){
					if(k == 1){
						rH = permutuj(k,t,h,phi);
						pi = 1;
						dopermutuj(pi,phi);
						//peknevypis(k,t,phi);
						phi->reset();
						nauto++;
						ncopf++;
					}
					else for(h = 0; h < n; h += k){
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
										else /*if(phi->rad()!=n/k)*/{
											
											//peknevypis(k,t,phi);
											ncopf++;
											//cout <<n<<" "<<k<<" "<<t<<" "<<h<<" "<<pi<<"\n";
										}
										//cekujorbity(k, phi);
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
	
	//printf("Total number of C.O.P.F. of C_%d is %d\n",n,ncopf);
	//printf("check subtotal of automorphisms of C_%d is %d\n",n,nauto);
	vysledok = ncopf;
	
	
	return vysledok;
}

bool coel(permutation *phi){
	bool ok = true;
	ok &= (phi->pi[1]!=1);
	for(int i=1; i<phi->cn; i++) ok &= (phi->pi[1] == (phi->pi[phi->cyklus[i]]));
	return ok;
}

bool coelorauto(permutation *phi){
	bool ok = true;
	//ok &= (phi->pi[1]!=1);
	for(int i=1; i<phi->cn; i++) ok &= (phi->pi[1] == (phi->pi[phi->cyklus[i]]));
	return ok;
}

permutation * umocni(permutation *phi, int k){
	int i, j, r;
	r = phi->rad();
	permutation * um = new permutation(phi->n);
	for(i = 0; i < phi->n; i++){
		um->predpis[i] = phi->predpis[i];
		um->pi[i] = phi->pi[i];
		for(j = 1; j < k; j++){
			um->pi[i] += phi->pi[um->predpis[i]];
			um->predpis[i] = phi->predpis[um->predpis[i]];
		}
		um->pi[i] = (um->pi[i]/k)%(r/k);
	}
	return um;
}

void sprav_zoznam(int n, int pocet, permutation ** zoznam, int copocet, permutation ** cozoznam, bool * copf){
	int i, j, k[pocet];
	int zabiich;
	permutation ** umocnene = new permutation *[pocet];
	for(i = 0; i < pocet; i++){
		k[i] = zoznam[i]->najdi_k();
		umocnene[i] = umocni(zoznam[i], k[i]);
	}
	for(i = 0; i < copocet; i++){
		zabiich = 0;
		for(j = 0; j < pocet; j++){
			if(porovnaj(umocnene[j], cozoznam[i])){
				zabiich++;
			}
		}
		if(zabiich){
			cout << "Z_" << n << "\n";
			cout << "\n--------------------------------------------------------------------\n";
			cozoznam[i]->vypis(false);
			if(copf[i]) cout << " copf\n";
			else cout << " auto\n";
			cozoznam[i]->vypispower(false);
			cout << " power functions\n\n";
			for(j = 0; j < pocet; j++){
				if(porovnaj(umocnene[j], cozoznam[i])){
					zoznam[j]->vypis(false);
					if(copf[i]) cout << " copf = skew ^ " << k[j] << "\n";
					else cout << " auto = skew ^ " << k[j] << "\n";
					zoznam[j]->vypispower(false);
					cout << " power functions\n\n";
				}
			}
		}
		if(zabiich) cout << "\n--------------------------------------------------------------------\n";
	}
	cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n";
}

int main(){
	ifstream f;
	f.open("zoznam/zoznam");
	int n, pocet, cyklov, dlzka, naozaj, conaozaj, aunaozaj, dok, dokco, dokau, a, b , z;
	permutation *zoznam[300];
	permutation *cozoznam[300];
	bool copf[300];
	permutation *phi;
	dok = dokco = dokau = 0;
    /*
	for(int i = 2; i<=52; i++){
		naozaj = 0;
		conaozaj = 0;
		aunaozaj = 0;
		f >> n >> pocet;
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
			if(!coelorauto(phi)){
				//zoznam[naozaj++] = phi;
				//phi = new permutation(n);
				naozaj++;
			}
			else{
				//copf[conaozaj] = coel(phi);
				//cozoznam[conaozaj++] = phi;
				//phi = new permutation(n);
				if(coel(phi)) conaozaj++;
				else aunaozaj++;
			}
		}
		
		//najdi(n, naozaj, zoznam);
		
		//sprav_zoznam(n, naozaj, zoznam, conaozaj, cozoznam, copf);
		
		dok+=naozaj;
		dokau+=aunaozaj;
		dokco+=conaozaj;
		
		printf("  Z_%2d auto %3d; copf %3d; zvysok %3d; dokopy %4d\n",n,aunaozaj,conaozaj,naozaj,naozaj+aunaozaj+conaozaj);
		
		//for(int j=0; j<naozaj; j++) delete zoznam[j];
		//for(int j=0; j<conaozaj; j++) delete cozoznam[j];
	}
	printf("dokopy auto %3d; copf %3d; zvysok %3d; dokopy %4d\n",dokau,dokco,dok,dok+dokau+dokco);
	*/
	f.close();

	int max = 0, nmax=1, sucet = 0;
	for(n=2;n<=500;n++){
		//cout << "n = " <<n<<"\n\n";
		pocet = najdi(n);
		sucet += pocet;
		if(max < pocet) {
			max = pocet;
			nmax = n;
		}
		cout << n << " " << pocet << "\n";
		//cout <<"\n..............................................................................\n\n";
		//if(pocet) cout << "Pre n="<<n<<" je pocet " << pocet<<"\n";
	}
	cout << "dokopy " << sucet << "\n";
	cout << "max: "<< nmax << " " << max << "\n";

}
