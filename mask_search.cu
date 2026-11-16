// mask_search.cu : fichier unique (core.h fusionné) pour compilation sur Kaggle.
// nvcc -O3 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_75,code=sm_75 -o mask_search mask_search.cu
// Usage: ./mask_search --mask 101x01x... --hash160 <hex40> [--device 0]
#include <cuda_runtime.h>
#include <chrono>
// core.h : GPU va CPU uchun umumiy kod (secp256k1 + SHA256 + RIPEMD160 + maska qidiruv "worker").
// Hech qanday tashqi kutubxona kerak emas. CPU testi aynan shu kodni ishlatadi.
#pragma once
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef __CUDACC__
#define HD __host__ __device__ __forceinline__
#else
#define HD inline
#endif
typedef uint32_t u32; typedef uint64_t u64;

#define LOWB 7                 // past x bitlar soni (batch = 2^LOWB = 128 nuqta)
#define MAXM (1<<LOWB)
#define MAXHI 64

// ---------------- 256-bit maydon elementi (little-endian 8x32) ----------------
struct fe { u32 d[8]; };
HD void fe_zero(fe&r){ for(int i=0;i<8;i++) r.d[i]=0; }
HD void fe_one(fe&r){ fe_zero(r); r.d[0]=1; }
HD void fe_setp(fe&r){ r.d[0]=0xFFFFFC2Fu; r.d[1]=0xFFFFFFFEu; for(int i=2;i<8;i++) r.d[i]=0xFFFFFFFFu; }
HD int fe_cmp(const fe&a,const fe&b){ for(int i=7;i>=0;i--){ if(a.d[i]!=b.d[i]) return a.d[i]>b.d[i]?1:-1; } return 0; }
HD int fe_is_zero(const fe&a){ u32 o=0; for(int i=0;i<8;i++) o|=a.d[i]; return o==0; }
HD u32 add256(fe&r,const fe&a,const fe&b){ u64 c=0; for(int i=0;i<8;i++){ c+=(u64)a.d[i]+b.d[i]; r.d[i]=(u32)c; c>>=32; } return (u32)c; }
HD u32 sub256(fe&r,const fe&a,const fe&b){ u64 br=0; for(int i=0;i<8;i++){ u64 t=(u64)a.d[i]-b.d[i]-br; r.d[i]=(u32)t; br=(t>>32)&1; } return (u32)br; }

HD void fe_add(fe&r,const fe&a,const fe&b){
  fe p; fe_setp(p);
  u32 c=add256(r,a,b);
  if(c){ fe k; fe_zero(k); k.d[0]=977; k.d[1]=1; add256(r,r,k); }   // 2^256 = 2^32+977 (mod p)
  else if(fe_cmp(r,p)>=0) sub256(r,r,p);
}
HD void fe_sub(fe&r,const fe&a,const fe&b){
  u32 br=sub256(r,a,b);
  if(br){ fe p; fe_setp(p); add256(r,r,p); }
}
HD void fe_neg(fe&r,const fe&a){ if(fe_is_zero(a)){ r=a; return; } fe p; fe_setp(p); sub256(r,p,a); }

HD void fe_mul(fe&r,const fe&a,const fe&b){
  u32 t[16];
  for(int i=0;i<16;i++) t[i]=0;
  for(int i=0;i<8;i++){
    u64 c=0;
    for(int j=0;j<8;j++){ u64 v=(u64)a.d[i]*b.d[j]+t[i+j]+c; t[i+j]=(u32)v; c=v>>32; }
    t[i+8]=(u32)c;
  }
  fe res; u64 c=0;
  for(int i=0;i<8;i++){
    u64 v=(u64)t[i]+(u64)t[8+i]*977u+c;
    if(i>0) v+=t[7+i];
    res.d[i]=(u32)v; c=v>>32;
  }
  u64 top=(u64)t[15]+c;                       // qiymat = res + top*2^256
  while(top){
    u64 m=top*977u; fe k; fe_zero(k);
    k.d[0]=(u32)m;
    u64 x1=(m>>32)+(u32)top;  k.d[1]=(u32)x1;
    u64 x2=(x1>>32)+(top>>32); k.d[2]=(u32)x2; k.d[3]=(u32)(x2>>32);
    top=add256(res,res,k);
  }
  fe p; fe_setp(p);
  if(fe_cmp(res,p)>=0) sub256(res,res,p);
  r=res;
}
HD void fe_sqr(fe&r,const fe&a){ fe_mul(r,a,a); }

// a^(p-2) mod p (Fermat). p-2 = FFFF...FFFE FFFFFC2D
HD void fe_inv(fe&r,const fe&a){
  fe res; fe_one(res);
  for(int i=255;i>=0;i--){
    fe_sqr(res,res);
    int k=i>>5; u32 e = (k==0)?0xFFFFFC2Du : (k==1)?0xFFFFFFFEu : 0xFFFFFFFFu;
    if((e>>(i&31))&1u) fe_mul(res,res,a);
  }
  r=res;
}

// ---------------- Egri chiziq (affin) ----------------
struct ge { fe x,y; int inf; };
HD void ge_set_inf(ge&r){ fe_zero(r.x); fe_zero(r.y); r.inf=1; }
HD void ge_add(ge&r,const ge&a,const ge&b){
  if(a.inf){ r=b; return; }
  if(b.inf){ r=a; return; }
  fe lam,t,u;
  if(fe_cmp(a.x,b.x)==0){
    fe s; fe_add(s,a.y,b.y);
    if(fe_is_zero(s)){ ge_set_inf(r); return; }
    fe_sqr(t,a.x); fe_add(u,t,t); fe_add(t,u,t);          // 3x^2
    fe_add(u,a.y,a.y); fe_inv(u,u); fe_mul(lam,t,u);
  } else {
    fe_sub(t,b.y,a.y); fe_sub(u,b.x,a.x); fe_inv(u,u); fe_mul(lam,t,u);
  }
  fe x3,y3;
  fe_sqr(x3,lam); fe_sub(x3,x3,a.x); fe_sub(x3,x3,b.x);
  fe_sub(t,a.x,x3); fe_mul(y3,lam,t); fe_sub(y3,y3,a.y);
  r.x=x3; r.y=y3; r.inf=0;
}
HD void ge_neg(ge&r,const ge&a){ r=a; if(!a.inf) fe_neg(r.y,a.y); }

// ---------------- SHA256 (33 baytlik siqilgan pubkey, bitta blok) ----------------
HD u32 ror32(u32 x,int n){ return (x>>n)|(x<<(32-n)); }
HD u32 rol32(u32 x,int n){ return (x<<n)|(x>>(32-n)); }
HD u32 bswap32(u32 x){ return (x>>24)|((x>>8)&0xff00u)|((x<<8)&0xff0000u)|(x<<24); }

HD void sha256_pub33(u32 out[8], u32 prefix, const fe&x){
  static const u32 K[64]={
  0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
  0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
  0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
  0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
  0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
  0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
  0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
  0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
  u32 w[64], xw[8];
  for(int k=0;k<8;k++) xw[k]=x.d[7-k];                 // X big-endian so'zlari
  w[0]=(prefix<<24)|(xw[0]>>8);
  for(int k=1;k<8;k++) w[k]=(xw[k-1]<<24)|(xw[k]>>8);
  w[8]=(xw[7]<<24)|0x00800000u;
  for(int k=9;k<15;k++) w[k]=0;
  w[15]=264;                                            // 33 bayt = 264 bit
  for(int i=16;i<64;i++){
    u32 s0=ror32(w[i-15],7)^ror32(w[i-15],18)^(w[i-15]>>3);
    u32 s1=ror32(w[i-2],17)^ror32(w[i-2],19)^(w[i-2]>>10);
    w[i]=w[i-16]+s0+w[i-7]+s1;
  }
  u32 a=0x6a09e667,b=0xbb67ae85,c=0x3c6ef372,d=0xa54ff53a,e=0x510e527f,f=0x9b05688c,g=0x1f83d9ab,h=0x5be0cd19;
  for(int i=0;i<64;i++){
    u32 S1=ror32(e,6)^ror32(e,11)^ror32(e,25);
    u32 ch=(e&f)^((~e)&g);
    u32 t1=h+S1+ch+K[i]+w[i];
    u32 S0=ror32(a,2)^ror32(a,13)^ror32(a,22);
    u32 mj=(a&b)^(a&c)^(b&c);
    u32 t2=S0+mj;
    h=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
  }
  out[0]=a+0x6a09e667; out[1]=b+0xbb67ae85; out[2]=c+0x3c6ef372; out[3]=d+0xa54ff53a;
  out[4]=e+0x510e527f; out[5]=f+0x9b05688c; out[6]=g+0x1f83d9ab; out[7]=h+0x5be0cd19;
}

// ---------------- RIPEMD160 (32 baytlik kirish, bitta blok) ----------------
HD u32 rmd_f(int j,u32 x,u32 y,u32 z){
  switch(j){ case 0: return x^y^z; case 1: return (x&y)|(~x&z); case 2: return (x|~y)^z;
             case 3: return (x&z)|(y&~z); default: return x^(y|~z); }
}
HD void ripemd160_32(u32 out[5], const u32 in_be[8]){   // in_be: SHA256 chiqishi (big-endian so'zlar)
  static const unsigned char R1[80]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15, 7,4,13,1,10,6,15,3,12,0,9,5,2,14,11,8,
    3,10,14,4,9,15,8,1,2,7,0,6,13,11,5,12, 1,9,11,10,0,8,12,4,13,3,7,15,14,5,6,2, 4,0,5,9,7,12,2,10,14,1,3,8,11,6,15,13};
  static const unsigned char R2[80]={5,14,7,0,9,2,11,4,13,6,15,8,1,10,3,12, 6,11,3,7,0,13,5,10,14,15,8,12,4,9,1,2,
    15,5,1,3,7,14,6,9,11,8,12,2,10,0,4,13, 8,6,4,1,3,11,15,0,5,12,2,13,9,7,10,14, 12,15,10,4,1,5,8,7,6,2,13,14,0,3,9,11};
  static const unsigned char S1[80]={11,14,15,12,5,8,7,9,11,13,14,15,6,7,9,8, 7,6,8,13,11,9,7,15,7,12,15,9,11,7,13,12,
    11,13,6,7,14,9,13,15,14,8,13,6,5,12,7,5, 11,12,14,15,14,15,9,8,9,14,5,6,8,6,5,12, 9,15,5,11,6,8,13,12,5,12,13,14,11,8,5,6};
  static const unsigned char S2[80]={8,9,9,11,13,15,15,5,7,7,8,11,14,14,12,6, 9,13,15,7,12,8,9,11,7,7,12,7,6,15,13,11,
    9,7,15,11,8,6,6,14,12,13,5,14,13,13,7,5, 15,5,8,11,14,14,6,14,6,9,12,9,12,5,15,8, 8,5,12,9,12,5,14,6,8,13,6,5,15,13,11,11};
  static const u32 K1[5]={0x00000000u,0x5A827999u,0x6ED9EBA1u,0x8F1BBCDCu,0xA953FD4Eu};
  static const u32 K2[5]={0x50A28BE6u,0x5C4DD124u,0x6D703EF3u,0x7A6D76E9u,0x00000000u};
  u32 X[16];
  for(int i=0;i<8;i++) X[i]=bswap32(in_be[i]);
  X[8]=0x80; for(int i=9;i<14;i++) X[i]=0; X[14]=256; X[15]=0;
  u32 h0=0x67452301u,h1=0xEFCDAB89u,h2=0x98BADCFEu,h3=0x10325476u,h4=0xC3D2E1F0u;
  u32 a=h0,b=h1,c=h2,d=h3,e=h4, a2=h0,b2=h1,c2=h2,d2=h3,e2=h4;
  for(int j=0;j<80;j++){
    int r=j>>4;
    u32 t=rol32(a+rmd_f(r,b,c,d)+X[R1[j]]+K1[r],S1[j])+e;
    a=e; e=d; d=rol32(c,10); c=b; b=t;
    t=rol32(a2+rmd_f(4-r,b2,c2,d2)+X[R2[j]]+K2[r],S2[j])+e2;
    a2=e2; e2=d2; d2=rol32(c2,10); c2=b2; b2=t;
  }
  u32 t=h1+c+d2; h1=h2+d+e2; h2=h3+e+a2; h3=h4+a+b2; h4=h0+b+c2; h0=t;
  out[0]=h0; out[1]=h1; out[2]=h2; out[3]=h3; out[4]=h4;   // digest baytlari = shu so'zlarning little-endian ko'rinishi
}
HD void hash160_pub(u32 out[5], u32 prefix, const fe&x){ u32 h[8]; sha256_pub33(h,prefix,x); ripemd160_32(out,h); }

// ---------------- Ish tavsifi va worker ----------------
struct Job {
  fe base; int L, nlo, nhi;
  int poslo[LOWB], poshi[MAXHI];
  ge B0;                       // base * G
  ge Q[MAXM];                  // Q[a] = (poslo bitlari a bo'yicha yig'indi) * G ; Q[0] = cheksiz
  ge Thi[MAXHI];               // 2^poshi[j] * G
  u32 target[5];               // hash160, 5 ta little-endian so'z
};
struct Result { int found; u64 s; u32 a; };

HD int ctz64(u64 v){ int n=0; while(!(v&1)){ v>>=1; n++; } return n; }

HD void check_hit(const Job*J, Result*R, const fe&x, const fe&y, u64 s, u32 a){
  u32 h[5]; hash160_pub(h, 2u+(y.d[0]&1u), x);
  if(h[0]==J->target[0]&&h[1]==J->target[1]&&h[2]==J->target[2]&&h[3]==J->target[3]&&h[4]==J->target[4]){
#ifdef __CUDA_ARCH__
    if(atomicCAS(&R->found,0,1)==0){ R->s=s; R->a=a; __threadfence(); }
#else
    if(!R->found){ R->found=1; R->s=s; R->a=a; }
#endif
  }
}

// Bitta B nuqta uchun 2^nlo ta (B + Q[a]) ni tekshirish: bitta inversiya (Montgomery batch)
HD void check_batch(const Job*J, Result*R, const ge&B, u64 s){
  const int M=1<<J->nlo;
  check_hit(J,R,B.x,B.y,s,0);                               // a = 0: B ning o'zi
  if(M==1) return;
  fe dx[MAXM], pr[MAXM], acc; fe_one(acc);
  for(int a=1;a<M;a++){
    fe d; fe_sub(d,J->Q[a].x,B.x);
    if(fe_is_zero(d)) fe_one(d);                            // (amalda bo'lmaydi), keyin o'tkazib yuboriladi
    dx[a]=d; fe_mul(acc,acc,d); pr[a]=acc;
  }
  fe inv; fe_inv(inv,pr[M-1]);
  for(int a=M-1;a>=1;a--){
    fe invd; if(a>1) fe_mul(invd,inv,pr[a-1]); else invd=inv;
    fe_mul(inv,inv,dx[a]);
    if(fe_cmp(J->Q[a].x,B.x)==0) continue;
    fe lam,t,x3,y3;
    fe_sub(t,J->Q[a].y,B.y); fe_mul(lam,t,invd);
    fe_sqr(x3,lam); fe_sub(x3,x3,B.x); fe_sub(x3,x3,J->Q[a].x);
    fe_sub(t,B.x,x3); fe_mul(y3,lam,t); fe_sub(y3,y3,B.y);
    check_hit(J,R,x3,y3,s,(u32)a);
  }
}

// Bitta thread/ish birligi: s0 dan boshlab 'per' ta yuqori-indeks qadam (Gray code)
HD void worker(const Job*J, Result*R, u64 s0, u64 per, u64 S){
  if(s0>=S) return;
  u64 s1=s0+per; if(s1>S) s1=S;
  ge B=J->B0; u64 g=s0^(s0>>1);
  for(int j=0;j<J->nhi;j++) if((g>>j)&1) ge_add(B,B,J->Thi[j]);
  for(u64 s=s0;s<s1;s++){
    if(*(volatile int*)&R->found) return;
    check_batch(J,R,B,s);
    if(s+1<s1){
      int j=ctz64(s+1);
      if((g>>j)&1){ ge n; ge_neg(n,J->Thi[j]); ge_add(B,B,n); } else ge_add(B,B,J->Thi[j]);
      g^=(u64)1<<j;
    }
  }
}

// ---------------- Host: maskadan Job yasash va natijani tekshirish ----------------
inline void fe_setbit(fe&a,int b){ a.d[b>>5]|=1u<<(b&31); }
inline int  fe_getbit(const fe&a,int b){ return (a.d[b>>5]>>(b&31))&1; }

inline int build_job(const char*mask,const char*hex40,Job&J,char*err){
  int L=(int)strlen(mask);
  if(L<2||L>250){ strcpy(err,"maska uzunligi noto'g'ri"); return 0; }
  if(mask[0]!='1'){ strcpy(err,"maskaning birinchi belgisi '1' bo'lishi shart"); return 0; }
  memset(&J,0,sizeof(J));
  J.L=L; fe_zero(J.base);
  int pos[256], nx=0;
  for(int i=L-1;i>=0;i--){                 // bit = L-1-i, pastdan yuqoriga
    int bit=L-1-i; char c=mask[i];
    if(c=='1') fe_setbit(J.base,bit);
    else if(c=='x'||c=='X') pos[nx++]=bit;
    else if(c!='0'){ strcpy(err,"maskada faqat 1, 0, x bo'lishi mumkin"); return 0; }
  }
  J.nlo = nx<LOWB?nx:LOWB; J.nhi = nx-J.nlo;
  if(J.nhi>MAXHI||J.nhi>62){ strcpy(err,"x lar juda ko'p"); return 0; }
  for(int i=0;i<J.nlo;i++) J.poslo[i]=pos[i];
  for(int j=0;j<J.nhi;j++) J.poshi[j]=pos[J.nlo+j];
  if(strlen(hex40)!=40){ strcpy(err,"hash160 40 ta hex belgi bo'lishi kerak"); return 0; }
  unsigned char hb[20];
  for(int i=0;i<20;i++){ unsigned v; if(sscanf(hex40+2*i,"%2x",&v)!=1){ strcpy(err,"hash160 hex xato"); return 0;} hb[i]=(unsigned char)v; }
  for(int i=0;i<5;i++) J.target[i]=hb[4*i]|(hb[4*i+1]<<8)|(hb[4*i+2]<<16)|((u32)hb[4*i+3]<<24);
  // D[k] = 2^k * G
  static ge D[256];
  const u32 gx[8]={0x16F81798,0x59F2815B,0x2DCE28D9,0x029BFCDB,0xCE870B07,0x55A06295,0xF9DCBBAC,0x79BE667E};
  const u32 gy[8]={0xFB10D4B8,0x9C47D08F,0xA6855419,0xFD17B448,0x0E1108A8,0x5DA4FBFC,0x26A3C465,0x483ADA77};
  for(int i=0;i<8;i++){ D[0].x.d[i]=gx[i]; D[0].y.d[i]=gy[i]; } D[0].inf=0;
  for(int k=1;k<L;k++) ge_add(D[k],D[k-1],D[k-1]);
  ge_set_inf(J.B0);
  for(int b=0;b<L;b++) if(fe_getbit(J.base,b)) ge_add(J.B0,J.B0,D[b]);
  ge_set_inf(J.Q[0]);
  for(int a=1;a<(1<<J.nlo);a++){ int i=ctz64(a); ge_add(J.Q[a],J.Q[a&(a-1)],D[J.poslo[i]]); }
  for(int j=0;j<J.nhi;j++) J.Thi[j]=D[J.poshi[j]];
  return 1;
}

inline void key_from(const Job&J,u64 s,u32 a,fe&key){
  key=J.base; u64 g=s^(s>>1);
  for(int j=0;j<J.nhi;j++) if((g>>j)&1) fe_setbit(key,J.poshi[j]);
  for(int i=0;i<J.nlo;i++) if((a>>i)&1) fe_setbit(key,J.poslo[i]);
}
// mustaqil tekshiruv: key*G -> hash160 == target ?
inline int verify_key(const Job&J,const fe&key){
  ge P; ge_set_inf(P);
  static ge D[256];
  const u32 gx[8]={0x16F81798,0x59F2815B,0x2DCE28D9,0x029BFCDB,0xCE870B07,0x55A06295,0xF9DCBBAC,0x79BE667E};
  const u32 gy[8]={0xFB10D4B8,0x9C47D08F,0xA6855419,0xFD17B448,0x0E1108A8,0x5DA4FBFC,0x26A3C465,0x483ADA77};
  for(int i=0;i<8;i++){ D[0].x.d[i]=gx[i]; D[0].y.d[i]=gy[i]; } D[0].inf=0;
  for(int k=1;k<J.L;k++) ge_add(D[k],D[k-1],D[k-1]);
  for(int b=0;b<J.L;b++) if(fe_getbit(key,b)) ge_add(P,P,D[b]);
  if(P.inf) return 0;
  u32 h[5]; hash160_pub(h,2u+(P.y.d[0]&1u),P.x);
  for(int i=0;i<5;i++) if(h[i]!=J.target[i]) return 0;
  return 1;
}
inline void key_hex(const fe&k,char*out){ for(int i=7;i>=0;i--){ sprintf(out,"%08x",k.d[i]); out+=8; } *out=0; }


__constant__ Job c_job;

__global__ void __launch_bounds__(128) kern(Result* R, u64 start, u64 per, u64 S){
  u64 tid = (u64)blockIdx.x*blockDim.x + threadIdx.x;
  worker(&c_job, R, start + tid*per, per, S);
}
#define CK(x) do{ cudaError_t e=(x); if(e!=cudaSuccess){ printf("CUDA xato %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); return 3; } }while(0)

int main(int argc,char**argv){
  const char *mask=0,*hash=0; int blocks=512, threads=128, dev=0; u64 per=64;
  for(int i=1;i<argc;i++){
    if(!strcmp(argv[i],"--mask")&&i+1<argc) mask=argv[++i];
    else if(!strcmp(argv[i],"--hash160")&&i+1<argc) hash=argv[++i];
    else if(!strcmp(argv[i],"--blocks")&&i+1<argc) blocks=atoi(argv[++i]);
    else if(!strcmp(argv[i],"--per")&&i+1<argc) per=strtoull(argv[++i],0,10);
    else if(!strcmp(argv[i],"--device")&&i+1<argc) dev=atoi(argv[++i]);
  }
  if(!mask||!hash){ printf("usage: %s --mask MASK --hash160 HEX40 [--blocks N] [--per N] [--device N]\n",argv[0]); return 2; }
  Job* J=(Job*)malloc(sizeof(Job)); char err[128];
  if(!build_job(mask,hash,*J,err)){ printf("XATO: %s\n",err); return 2; }
  u64 S=(u64)1<<J->nhi, keys=(u64)1<<(J->nhi+J->nlo);
  printf("maska uzunligi=%d, x=%d, jami kalit=2^%d, yuqori qadamlar=2^%d\n",J->L,J->nhi+J->nlo,J->nhi+J->nlo,J->nhi);
  CK(cudaSetDevice(dev));
  CK(cudaMemcpyToSymbol(c_job,J,sizeof(Job)));
  Result h={0,0,0}; Result* dR; CK(cudaMalloc(&dR,sizeof(Result))); CK(cudaMemcpy(dR,&h,sizeof(h),cudaMemcpyHostToDevice));
  u64 T=(u64)blocks*threads, chunk=T*per, done=0;
  auto t0=std::chrono::steady_clock::now();
  for(u64 start=0; start<S; start+=chunk){
    kern<<<blocks,threads>>>(dR,start,per,S);
    CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
    CK(cudaMemcpy(&h,dR,sizeof(h),cudaMemcpyDeviceToHost));
    done = (start+chunk<S?start+chunk:S);
    double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    printf("\r%.4f%%  %.1f Mkey/s   ", 100.0*done/S, (double)done*(1<<J->nlo)/sec/1e6); fflush(stdout);
    if(h.found) break;
  }
  printf("\n");
  if(!h.found){ printf("NOTFOUND: bu maskada kalit yo'q (pattern noto'g'ri yoki hash mos emas)\n"); return 1; }
  fe key; key_from(*J,h.s,h.a,key); char hx[80]; key_hex(key,hx);
  int ok=verify_key(*J,key);                    // CPU'da mustaqil qayta tekshirish
  printf("FOUND! private key (hex) = %s\nCPU tekshiruvi: %s\n",hx, ok?"TASDIQLANDI":"TASDIQLANMADI (xato!)");
  FILE*f=fopen("FOUND.txt","a"); if(f){ fprintf(f,"%s\n",hx); fclose(f); }
  return ok?0:4;
}
