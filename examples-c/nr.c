#include <math.h>
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(float data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

void fourn(float data[], unsigned long nn[], int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}

#undef SWAP

void realft(float data[], unsigned long n, int isign)
{
	void four1(float data[], unsigned long nn, int isign);
	unsigned long i,i1,i2,i3,i4,np3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(data,n>>1,-1);
	}
}



void rlft3(float ***data, float **speq, unsigned long nn1, unsigned long nn2,
	unsigned long nn3, int isign)
{
	void nrerror(char error_text[]);
	unsigned long i1,i2,i3,j1,j2,j3,nn[4],ii3;
	double theta,wi,wpi,wpr,wr,wtemp;
	float c1,c2,h1r,h1i,h2r,h2i;

	if (1+&data[nn1][nn2][nn3]-&data[1][1][1] != nn1*nn2*nn3)
		nrerror("rlft3: problem with dimensions or contiguity of data array\n");
	c1=0.5;
	c2 = -0.5*isign;
	theta=isign*(6.28318530717959/nn3);
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	nn[1]=nn1;
	nn[2]=nn2;
	nn[3]=nn3 >> 1;
	if (isign == 1) {
		fourn(&data[1][1][1]-1,nn,3,isign);
		for (i1=1;i1<=nn1;i1++)
			for (i2=1,j2=0;i2<=nn2;i2++) {
				speq[i1][++j2]=data[i1][i2][1];
				speq[i1][++j2]=data[i1][i2][2];
			}
	}
	for (i1=1;i1<=nn1;i1++) {
		j1=(i1 != 1 ? nn1-i1+2 : 1);
		wr=1.0;
		wi=0.0;
		for (ii3=1,i3=1;i3<=(nn3>>2)+1;i3++,ii3+=2) {
			for (i2=1;i2<=nn2;i2++) {
				if (i3 == 1) {
					j2=(i2 != 1 ? ((nn2-i2)<<1)+3 : 1);
					h1r=c1*(data[i1][i2][1]+speq[j1][j2]);
					h1i=c1*(data[i1][i2][2]-speq[j1][j2+1]);
					h2i=c2*(data[i1][i2][1]-speq[j1][j2]);
					h2r= -c2*(data[i1][i2][2]+speq[j1][j2+1]);
					data[i1][i2][1]=h1r+h2r;
					data[i1][i2][2]=h1i+h2i;
					speq[j1][j2]=h1r-h2r;
					speq[j1][j2+1]=h2i-h1i;
				} else {
					j2=(i2 != 1 ? nn2-i2+2 : 1);
					j3=nn3+3-(i3<<1);
					h1r=c1*(data[i1][i2][ii3]+data[j1][j2][j3]);
					h1i=c1*(data[i1][i2][ii3+1]-data[j1][j2][j3+1]);
					h2i=c2*(data[i1][i2][ii3]-data[j1][j2][j3]);
					h2r= -c2*(data[i1][i2][ii3+1]+data[j1][j2][j3+1]);
					data[i1][i2][ii3]=h1r+wr*h2r-wi*h2i;
					data[i1][i2][ii3+1]=h1i+wr*h2i+wi*h2r;
					data[j1][j2][j3]=h1r-wr*h2r+wi*h2i;
					data[j1][j2][j3+1]= -h1i+wr*h2i+wi*h2r;
				}
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
	}
	if (isign == -1)
		fourn(&data[1][1][1]-1,nn,3,isign);
}

void sinft(float y[], int n)
{
	void realft(float data[], unsigned long n, int isign);
	int j,n2=n+2;
	float sum,y1,y2;
	double theta,wi=0.0,wr=1.0,wpi,wpr,wtemp;

	theta=3.14159265358979/(double) n;
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	y[1]=0.0;
	for (j=2;j<=(n>>1)+1;j++) {
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
		y1=wi*(y[j]+y[n2-j]);
		y2=0.5*(y[j]-y[n2-j]);
		y[j]=y1+y2;
		y[n2-j]=y1-y2;
	}
	realft(y,n,1);
	y[1]*=0.5;
	sum=y[2]=0.0;
	for (j=1;j<=n-1;j+=2) {
		sum += y[j];
		y[j]=y[j+1];
		y[j+1]=sum;
	}
}

#define PI 3.141592653589793

void cosft1(float y[], int n)
{
	void realft(float data[], unsigned long n, int isign);
	int j,n2;
	float sum,y1,y2;
	double theta,wi=0.0,wpi,wpr,wr=1.0,wtemp;

	theta=PI/n;
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	sum=0.5*(y[1]-y[n+1]);
	y[1]=0.5*(y[1]+y[n+1]);
	n2=n+2;
	for (j=2;j<=(n>>1);j++) {
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
		y1=0.5*(y[j]+y[n2-j]);
		y2=(y[j]-y[n2-j]);
		y[j]=y1-wi*y2;
		y[n2-j]=y1+wi*y2;
		sum += wr*y2;
	}
	realft(y,n,1);
	y[n+1]=y[2];
	y[2]=sum;
	for (j=4;j<=n;j+=2) {
		sum += y[j];
		y[j]=sum;
	}
}
#undef PI

#define PI 3.141592653589793

void cosft2(float y[], int n, int isign)
{
	void realft(float data[], unsigned long n, int isign);
	int i;
	float sum,sum1,y1,y2,ytemp;
	double theta,wi=0.0,wi1,wpi,wpr,wr=1.0,wr1,wtemp;

	theta=0.5*PI/n;
	wr1=cos(theta);
	wi1=sin(theta);
	wpr = -2.0*wi1*wi1;
	wpi=sin(2.0*theta);
	if (isign == 1) {
		for (i=1;i<=n/2;i++) {
			y1=0.5*(y[i]+y[n-i+1]);
			y2=wi1*(y[i]-y[n-i+1]);
			y[i]=y1+y2;
			y[n-i+1]=y1-y2;
			wr1=(wtemp=wr1)*wpr-wi1*wpi+wr1;
			wi1=wi1*wpr+wtemp*wpi+wi1;
		}
		realft(y,n,1);
		for (i=3;i<=n;i+=2) {
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
			y1=y[i]*wr-y[i+1]*wi;
			y2=y[i+1]*wr+y[i]*wi;
			y[i]=y1;
			y[i+1]=y2;
		}
		sum=0.5*y[2];
		for (i=n;i>=2;i-=2) {
			sum1=sum;
			sum += y[i];
			y[i]=sum1;
		}
	} else if (isign == -1) {
		ytemp=y[n];
		for (i=n;i>=4;i-=2) y[i]=y[i-2]-y[i];
		y[2]=2.0*ytemp;
		for (i=3;i<=n;i+=2) {
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
			y1=y[i]*wr+y[i+1]*wi;
			y2=y[i+1]*wr-y[i]*wi;
			y[i]=y1;
			y[i+1]=y2;
		}
		realft(y,n,-1);
		for (i=1;i<=n/2;i++) {
			y1=y[i]+y[n-i+1];
			y2=(0.5/wi1)*(y[i]-y[n-i+1]);
			y[i]=0.5*(y1+y2);
			y[n-i+1]=0.5*(y1-y2);
			wr1=(wtemp=wr1)*wpr-wi1*wpi+wr1;
			wi1=wi1*wpr+wtemp*wpi+wi1;
		}
	}
}
#undef PI
