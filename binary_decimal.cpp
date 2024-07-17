/*
*/
#include "binary_decimal.h"
#include <assert.h>

int Act_Inv_Translation_assuming_PBC(int dec_state, int l1_, int l2_ ){

    int l2=l2_;
    int l1=l1_;

    //PBC b/w l2 position and l1 position is assumed
    /*state_1
    ( _ _ _ _ _ _ _ _   _ ) 0  1  1  0 _ _ _ _  0  1 ( _  _ _ _ _ _ _ _ _ _ _ _ )
      0 1 2           l1-1  l1 l1+1                l2 l2+1
    */
    /* To ======>
    /*state_2
    ( _ _ _ SAME   _ _ _   ) 1  1  0  0 _ _ _ _  1  0 ( _ _ _ _ _ _ SAME _ _ _ )
                             l1 l1+1                l2
    */


    int dec_state_new;

    dec_state_new = dec_state;

    for(int i=l1+1;i<=l2;i++){
        dec_state_new -= bit_value(dec_state,i)*(int)(pow(2,i-1)+0.5);
    }

    dec_state_new += (bit_value(dec_state,l1)*(int)(pow(2,l2)+0.5))
            - (bit_value(dec_state,l1)*(int)(pow(2,l1)+0.5));

    return dec_state_new;

}

int Act_Translation_assuming_PBC(int dec_state, int l1_, int l2_ ){

    int l2=l2_;
    int l1=l1_;

    //PBC b/w l2 position and l1 position is assumed
    /*state_1
    ( _ _ _ _ _ _ _ _   _ ) 0  1  1  0 _ _ _ _  0  1 ( _  _ _ _ _ _ _ _ _ _ _ _ )
      0 1 2           l1-1  l1 l1+1                l2 l2+1
    */
    /* To ======>
    /*state_2
    ( _ _ _ SAME   _ _ _   ) 1  0  1  1 _ _ _ _  _  0  ( _ _ _ _ _ _ SAME _ _ _ )
                             l1 l1+1                l2
    */
    int dec_state_new;

    dec_state_new = dec_state;

    for(int i=l1;i<l2;i++){
        dec_state_new += bit_value(dec_state,i)*(int)(pow(2,i)+0.5);
    }
    dec_state_new += (bit_value(dec_state,l2)*(int)(pow(2,l1)+0.5))
            - (bit_value(dec_state,l2)*(int)(pow(2,l2)+0.5));

    return dec_state_new;

}




int Act_Translation_2D_alongY_assuming_PBC(int dec_state, int Lx, int Ly, int ix){


    //lattice index l = ix + iy*Lx
    int l;

    int dec_state_new;
    dec_state_new = dec_state;

    for(int iy=0;iy<Ly-1;iy++){
        l=ix+(iy*Lx);
        dec_state_new += bit_value(dec_state,l)*(int)(pow(2,l)+0.5)*
                         (int)(pow(2,Lx) - 1 + 0.5);
    }
    l=ix+((Ly-1)*Lx);
    dec_state_new += (bit_value(dec_state,l)*(int)(pow(2,ix)+0.5))
            - (bit_value(dec_state,l)*(int)(pow(2,l)+0.5));

    return dec_state_new;

}


int Act_Translation_2D_alongX_assuming_PBC(int dec_state, int Lx, int Ly, int iy){


    //lattice index l = ix + iy*Lx
    int l1, l2;
    l1 = iy*Lx;
    l2 = (Lx-1) + iy*Lx;

    int dec_state_new;

    dec_state_new=Act_Translation_assuming_PBC(dec_state, l1, l2 );

    return dec_state_new;
}


void print_binary_of_decimal(int n){

    unsigned i,j;
    j=1<<31;
    assert (n < j);
    for (i = 1 << 31; i > 0; i = i / 2){
        if((n & i) !=0){
            cout<<"1";
        }
        else{
            cout<<"0";
        }
    }
    cout<<endl;
}

Mat_1_int decimal_to_binary(int n){

    Mat_1_int Dec;
    Dec.resize(32);
    unsigned i,j;
    j=1<<31;
    assert (n < j);
    int pos=31;
    for (i = 1 << 31; i > 0; i = i / 2){
        if((n & i) !=0){
            Dec[pos]=1;
        }
        else{
            Dec[pos]=0;
        }
        pos=pos-1;
    }


    return Dec;
}

int one_bits_in_bw(int n, int m, int Decimal){

    assert((n<32) && (m<32));

    int l ,lp;
    int n_ones, D_overlap, D_;
    int overlap;

    if(n>=m){
        l=n;
        lp=m;
    }
    else{
        l=m;
        lp=n;
    }

    if ( (l==lp) || (l==lp-1)){
        overlap=0;
    }
    else{
        n_ones=l-lp-1;
        D_ = (int) (pow(2,n_ones)-1+0.5);
        D_overlap = D_ << (lp+1);
        overlap = countCommonBits(D_overlap,Decimal);
    }

    return overlap;
}

int countCommonBits(int a,int b) {
    int n = 0;
    for (unsigned v = (unsigned)(a & b); v; v >>= 1) {
        n += 1 & v;
    }
    return n;
}
