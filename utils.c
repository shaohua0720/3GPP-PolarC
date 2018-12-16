#include "utils.h"

uint32_t DS1CA_polar_encoder(uint8_t *a,  //info bits of length A
			     uint32_t A,
			     uint8_t *crc_polynomial, //of length P+1
			     uint32_t P,
			     uint8_t *crc_scrambling, //RNTI of length 16
			     uint8_t *crc_interleaver, //of length K
			     uint32_t K,
			     uint32_t *info_bit, // of length N
			     uint32_t N,
			     uint32_t *rate_matching, //of length E
			     uint32_t E,
			     uint32_t *e //output of length E
			     )
{
  uint32_t i,j;
  uint8_t **G_P=(uint8_t **)malloc((A+P)*sizeof(uint8_t*));
  for(i=0;i<A+P;i++) {
    G_P[i]=(uint8_t *)malloc(P*sizeof(uint8_t));
  }
  get_crc_generator_matrix(A+P,crc_polynomial,P,G_P);

  // generate crc bits
  uint8_t *crc_bits=(uint8_t *)malloc((P)*sizeof(uint8_t));
  uint8_t *tmp_value=(uint8_t *)malloc((A+P)*sizeof(uint8_t));
  for(i=0;i<A+P;i++) {
    if (i<P)
      tmp_value[i] = 1;
    else
      tmp_value[i] = a[i-P];
  }

  for(i=0;i<P;i++) {
    crc_bits[i]=0;
    for(j=0;j<A+P;j++) {
      crc_bits[i] += tmp_value[j]*G_P[j][i];
    }
    crc_bits[i]=crc_bits[i]%2;
  }
  //free the mem
  for(i=0;i<A+P;i++) {
    free(G_P[i]);G_P[i]=NULL;
  }
  free(G_P);G_P=NULL;
  free(tmp_value);tmp_value=NULL;

  // scramble the RNTI
  for(i=P-16;i<P;i++) {
    crc_bits[i]=uint8_t_xor(crc_bits[i],crc_scrambling[i-P+16]);
  }

  // append the crc bits
  uint8_t *b=(uint8_t *)malloc((A+P)*sizeof(uint8_t));
  for (i=0;i<A+P;i++) {
    if(i<A)
      b[i]=a[i];
    else
      b[i]=crc_bits[i-A];
  }
  free(crc_bits);crc_bits=NULL;

  uint8_t *c = (uint8_t *)malloc((A+P)*sizeof(uint8_t));
  for(i=0;i<A+P;i++) {
    c[i] = b[crc_interleaver[i]-1];
  }
  free(b);b=NULL;
  uint32_t k=0;
  // place the info bits and crc bits.
  uint8_t *u=(uint8_t *)malloc(N*sizeof(uint8_t));
  for(i=0;i<N;i++) {
    u[i]=0;
  }
  for(i=0;i<N;i++) {
    if(info_bit[i]!=0) {
      u[i]=c[k];
      k++;
    }
  }
  free(c);c=NULL;

  // polar encoder
  uint32_t **G_N=(uint32_t **)malloc(N*sizeof(uint32_t*));
  for(i=0;i<N;i++) {
    G_N[i]=(uint32_t*)malloc(N*sizeof(uint32_t));
  }

  get_G_N(G_N,round(log2(N)));
  
#ifdef debug
  printf("\nvalue of G_N:\n");
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      printf("%d", G_N[i][j]);
    }
    printf("\n");
  }
#endif
  
  uint8_t *d= (uint8_t*)malloc(N*sizeof(N));
  for(i=0;i<N;i++) {
    d[i]=0;
  }
  
  for(i=0;i<N;i++) {
    d[i]=0;
    for (j=0;j<N;j++) {
      d[i]+=u[j]*G_N[j][i];
    }
    d[i]=d[i]%2;
  }

  for(i=0;i<N;i++) {
    free(G_N[i]);G_N[i]=NULL;
  }
  free(G_N);G_N=NULL;
  free(u);u=NULL;

  // rate match the data d
  for(i=0;i<E;i++) {
    e[i] = d[rate_matching[i]];
  }
  free(d); d=NULL;
  return 0;
}
// store the elements in A but not in B into C
uint8_t find_in_array(uint32_t *C, uint32_t N, uint32_t val)
{
  uint32_t i;
  for(i=0;i<N;i++) {
    if (val==C[i]) {
      return 1;
    }
  }
  return 0;
}
uint8_t setdiff(uint32_t *A,uint32_t a_l,uint32_t *B,uint32_t b_l, uint32_t **C, uint32_t *N)
{
  uint32_t i,j,k=0;
  uint32_t *C_t=(uint32_t *) malloc(a_l*sizeof(uint32_t));

  for(i=0;i<a_l;i++) {
    if (!find_in_array(B,b_l,A[i])) {
      C_t[k]=A[i];
      k++;
    }
  }
  if(k==0) {
    free(C_t);C_t=NULL;
    *N=0;*C=NULL;
  } else {
    *C = (uint32_t *)realloc(C_t,k*sizeof(uint32_t));
    *N = k;    
  }
  return 0;
}
uint8_t get_3GPP_info_bit_pattern(uint32_t I, uint32_t *Q_N, uint32_t N, uint32_t *rate_matching, uint32_t E, mode mode_t, uint32_t *info_bit)
{
  uint32_t n = round(log2(N));
  if (I>N)
    {
      printf("Unsupported block length!I<=N\n");
      return 1;
    }

  if(I>E)
    {
      printf("Unspported block length,I<=E!\n");
    }

  if(mode_t==repetition){
      if(E<N){
	printf("mode is not compatible with E\n");
	return 1;
      }
    }
  else if(mode_t==puncturing) {
    if (E>=N) {
      printf("mode is not compatible with E\n");
      return 1;
    }
  }
  else if (mode_t == shortening) {
    if (E>=N) {
      printf("mode is not compatible with E\n");
      return 1;
    }
  }
  else
    {
      printf("Unsupported mode.\n");
      return 1;
    }

  uint32_t *Q_tmp_N,*Q_tmp_N_p, Q_tmp_N_t;
  uint32_t i,j;
  uint32_t *array1toN=(uint32_t *)malloc(N*sizeof(uint32_t));
  for(i=0;i<N;i++)
    array1toN[i]=i;
  //compute the diff set
  setdiff(array1toN,N,rate_matching,E,&Q_tmp_N, &Q_tmp_N_t);
  free(array1toN);array1toN=NULL;
      
  if(mode_t == puncturing) {
    uint32_t append=0;
    if (E>= 3.0*N/4) {
      append=ceil(3.0*N/4-E/2.0);
      //realloc the mem for Q_tmp_N
    } else {
      append = ceil(9.0*N/16-E/4.0);
    }
    printf("\n append:%d\n",append);
    Q_tmp_N_p = (uint32_t *)realloc(Q_tmp_N,(Q_tmp_N_t+append)*sizeof(uint32_t));
    if (!Q_tmp_N) {
      printf("error!\n");
      return 1;
      
    }
      for(i=Q_tmp_N_t;i<Q_tmp_N_t+append;i++){
	Q_tmp_N_p[i]=i-Q_tmp_N_t;
      }
      Q_tmp_N_t=Q_tmp_N_t+append;
  }

  Q_tmp_N_p=Q_tmp_N;
  
    uint32_t *Q_Itmp_N=NULL, Q_Itmp_N_t=0;
  setdiff(Q_N,N,Q_tmp_N_p,Q_tmp_N_t,&Q_Itmp_N,&Q_Itmp_N_t);
  
  if (Q_tmp_N_p !=NULL) {
    free(Q_tmp_N_p);Q_tmp_N_p=NULL;
  }

  if(Q_Itmp_N_t < I) {
    printf("Unsupported block length: too many frozen bits\n");
  }

  uint32_t *Q_I_N =Q_Itmp_N+Q_Itmp_N_t-I; //Q_I_N of length I
  for(i=0;i<N;i++) {
    info_bit[i]=0;
  }
  for(i=0;i<I;i++) {
    info_bit[Q_I_N[i]]=1;
  }
  Q_I_N=NULL;
  free(Q_Itmp_N);Q_Itmp_N=NULL;
  
  return 0;
}
uint32_t get_3GPP_sequence_pattern(uint32_t *Q_N, uint32_t N)
{
  const uint16_t Q_Nmax[1024]=Q_Nmax_t;
  uint32_t N_t=round(log2(N));
  if(p_pow(2,N_t)!=N)
    {
      printf("N should be power of 2.\n");
      return 1;
    }
  uint32_t i,j=0;
  for(i=0;i<1024;i++)
    {
      if(Q_Nmax[i]<N)
	{
	  Q_N[j]=Q_Nmax[i];
	  j++;
	}
    }
  return 0;
}

uint32_t get_3GPP_rate_matching_pattern(uint32_t K,uint32_t N, uint32_t E, uint32_t *rate_matching_pattern, mode * mode_t)
{
  uint32_t i,j;
  uint32_t n=round(log2(N));
  printf("n:%d\n",n);
  if (n<5)
    {
      printf("polar_3gpp_c: unsupported block length, N should larger than 5.n");
      return 1;
    }
  uint32_t P[]={0,1,2,4,3,5,6,7,8,16,9,17,10,18,11,19,12,20,13,21,14,22,15,23,24,25,26,28,27,29,30,31};
  uint32_t *d=(uint32_t *)malloc(N*sizeof(uint32_t));
  for(i=0;i<N;i++)
    {
      d[i]=i+1;
      //      d[i]=i;
    }

  
  uint32_t *J=(uint32_t *)malloc(N*sizeof(uint32_t));
  uint32_t *y=(uint32_t *)malloc(N*sizeof(uint32_t));

  for(i=0;i<N;i++)
    {
      uint32_t idx=floor(32*i/N);
      J[i]=P[idx]*round(N/32)+(i%(N/32));
      y[i]=d[J[i]];
    }
  free(d);d=NULL;

  if (E>=N)
    {
      for(j=0;j<E;j++)
	{
	  rate_matching_pattern[j]=y[j%N]-1;
	}
      *mode_t=repetition;
    }
  else
    {
      if (1.0*K/E<=7.0/16)
       {
	for(j=0;j<E;j++) {
	    rate_matching_pattern[j]=y[j+N-E]-1;
	  }
	*mode_t=puncturing;
      } else {
      for(j=0;j<E;j++) {
	rate_matching_pattern[j]=y[j]-1;
      }
      *mode_t=shortening;
    }
  }
  free(J);J=NULL;
  free(y);y=NULL;
  return 0;
}

double log2(double n)
{
  return log(n)/log(2);
}

uint32_t get_3GPP_N(uint32_t K, uint32_t E,uint8_t n_max)
{
  uint32_t n1,n2,n;
  uint32_t val1= 9*p_pow(2,ceil(log2(E))-1)/8;
  if( (E<=val1) && (1.0*K/E <9.0/16))
    {
      n1=ceil(log2(E))-1;
    }else
    {
      n1=ceil(log2(E));
    }
  uint8_t n_min=5;
  n2=ceil(log2(K*8));
  n=max(n_min,min(min(n1,n2),n_max));
  return p_pow(2,n);
}

uint32_t p_pow(uint32_t base,uint32_t N)
{
  if(base==2)
      return 0x01<<N;
  uint32_t rs=1;
  while(N--)
    rs=rs*base;
  return rs;
}

uint8_t get_G_N(uint32_t **tmp,uint8_t N)
{
  //  printf("test p_pow:%d\n",p_p_pow(2,10));
  int i,j;
  // initialize a
  uint32_t **a=(uint32_t **)malloc(sizeof(uint32_t *)*2);
  a[0]=(uint32_t*)malloc(sizeof(uint32_t)*2);
  a[1]=(uint32_t*)malloc(sizeof(uint32_t)*2);
  a[0][0]=1;
  a[1][0]=1;
  a[0][1]=0;
  a[1][1]=1;
  
  uint32_t order=p_pow(2,N);

  tmp[0][0]=1;
  tmp[0][1]=0;
  tmp[1][0]=1;
  tmp[1][1]=1;

  uint32_t *C;
  
  for(i=2;i<=N;i++)
    {
      //alloc C
      uint32_t **C=(uint32_t **)malloc(sizeof(uint32_t*)*p_pow(2,i));
      for(j=0;j<p_pow(2,i);j++)
	{
	  C[j]=(uint32_t *)malloc(sizeof(uint32_t)*p_pow(2,i));
	}
      kronecker(a,2,2,tmp,p_pow(2,i-1),p_pow(2,i-1),C);
      //copy C to tmp
      uint32_t p,q;
      for(p=0;p<p_pow(2,i);p++)
	{
	  for(q=0;q<p_pow(2,i);q++)
	    tmp[p][q]=C[p][q];
	}
      //free C
      for(p=0;p<p_pow(2,i);p++)
	{
	  free(C[p]);C[p]=NULL;
	}
      free(C);C=NULL; 
    }

  //free a
  for(i=0;i<2;i++)
    {
      free(a[i]);a[i]=NULL;
    }
  free(a);a=NULL;
  
  #ifdef DEBUG
  for(i=0;i<order;i++)
    {
      for(j=0;j<order;j++)
	{
	  //	  printf("%d  ",tmp[i][j]);
	}
      //      printf("\n");
    }
  #endif
  return 0;
}

uint8_t kronecker(uint32_t **A, uint8_t ai, uint8_t aj,uint32_t **B, uint8_t bi, uint8_t bj, uint32_t **C)
{
  int i,j;
  for(i=0;i<ai*bi;i++)
    {
      for(j=0;j<aj*bj;j++)
      {
	uint8_t rowidx0=i/bi;
	uint8_t colidx0=j/bj;
	uint32_t value0=A[rowidx0][colidx0];
	uint8_t rowidx1=i%bi;
	uint8_t colidx1=j%bj;
	uint32_t value1=B[rowidx1][colidx1];
	C[i][j]=value0*value1;
       #ifdef DEBUG
	//printf("A[%d][%d]=%d  ",rowidx0,colidx0,value0);
	//printf("B[%d][%d]=%d  ",rowidx1,colidx1,value1);
	//printf("%d  ",C[i][j]);
       #endif
      }
  #ifdef DEBUG
      //printf("\n");
  #endif
    }
  return 0;
}
uint8_t get_3GPP_crc_interleaver_pattern(uint8_t K, uint8_t *Pi)
{
  uint8_t IL[224]=IL_kernel;
  int i,j=0;
  for(i=0;i<K;i++)
    {
      Pi[i]=0;
    }
  
  for(i=0;i<224;i++)
    {
      //      printf("%d ",IL[i]);
      if(IL[i+1]>=(224-K))
	{
	  Pi[j]=IL[i+1]-(224-K)+1;
	  j++;
	}
    }
  return 0;
}

uint8_t uint8_t_xor(uint8_t a, uint8_t b)
{
  return a^b;
}

uint32_t get_crc_generator_matrix(uint32_t size,
				  uint8_t *crc_init,
				  uint32_t p,
				  uint8_t **G_P // of size x p
				  )
{

  //  uint8_t  crc_init[25]={1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1};
  uint32_t i;
  int j;
  //  printf("ttttt;%d %d\n",size,p);
  // Fill the G_P
  for(i=0;i<p;i++)
    {
      G_P[size-1][i]=crc_init[i+1];
    }
  for(j=size-2;j>=0;--j)
    {
      if(j<0)
	break;
      //  printf("j value:%d\n",j);
      uint8_t *tmp_GP1=malloc(sizeof(uint8_t)*p);
      if(!tmp_GP1)
	{
	  printf("mem error!");
	  return 1;
	}
      for(i=0;i<p-1;i++)
	{
	  tmp_GP1[i]=G_P[j+1][i+1];
	}
      tmp_GP1[p-1]=0;

      uint8_t tmp_value=G_P[j+1][0];
      uint8_t *tmp_GP2=malloc(sizeof(uint8_t)*p);
      if(!tmp_GP2)
	{
	  printf("Mem error!\n");
	  return 1;
	}
      for(i=0;i<p;i++)
	{
	  tmp_GP2[i]=tmp_value*crc_init[i+1];
	}

      for(i=0;i<p;i++)
	{
	  G_P[j][i]=uint8_t_xor(tmp_GP1[i],tmp_GP2[i]);
	}
      //      free(tmp_GP1);free(tmp_GP2);
      //      tmp_GP1=NULL; tmp_GP2=NULL;
   }

    
  // for test
 #ifdef DEBUG
  //  printf("ffffffff:%d %d\n",size,p);
  for(i=0;i<size;i++)
    {
      for(j=0;j<p;j++)
	{
	  //	  printf("%d ",G_P[i][j]);
	}
      //          printf("\n");
    }
  #endif
  return 0;
}

uint8_t CRC_cal(uint8_t *crc,uint8_t *a, uint8_t size,uint8_t **G_P, uint8_t P)
{
  int i,j;
  // Check the input data
  if (size==0)
    return 0;
  for(i=0;i<size;i++)
    {
      if (a[i]!=0&a[i]!=1)
	{
	  printf("CRC input should be 0 or 1!\n");
	  return 1;
	}
    }
  // construct [ones(1,p) a]
  uint8_t *src=(uint8_t *)malloc(sizeof(uint8_t)*(size+P));
  
  for(i=0;i<size+P;i++)
    {
      if(i<P)
	src[i]=1;
      else
	{
	  src[i]=a[i-P];
	}
    }

  //vector-matrix multiplication
  for(i=0;i<P;i++)
    {
      int colum_sum=0;
      for(j=0;j<size+P;j++)
	{
	  colum_sum+=src[j]*G_P[j][i];
	}
      crc[i]=colum_sum%2;
    }
  free(src);
  src=NULL;
  return 0;
}
