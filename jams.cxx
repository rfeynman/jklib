void jam555(double * out,double * in1,double * in2)
{
      double tem[25];
      tem[ 0]=in1[ 0]*in2[ 0] + in1[ 5]*in2[ 1] + in1[10]*in2[ 2] + in1[15]*in2[ 3] + in1[20]*in2[ 4];
      tem[ 5]=in1[ 0]*in2[ 5] + in1[ 5]*in2[ 6] + in1[10]*in2[ 7] + in1[15]*in2[ 8] + in1[20]*in2[ 9];
      tem[10]=in1[ 0]*in2[10] + in1[ 5]*in2[11] + in1[10]*in2[12] + in1[15]*in2[13] + in1[20]*in2[14];
      tem[15]=in1[ 0]*in2[15] + in1[ 5]*in2[16] + in1[10]*in2[17] + in1[15]*in2[18] + in1[20]*in2[19];
      tem[20]=in1[ 0]*in2[20] + in1[ 5]*in2[21] + in1[10]*in2[22] + in1[15]*in2[23] + in1[20]*in2[24];
      tem[ 1]=in1[ 1]*in2[ 0] + in1[ 6]*in2[ 1] + in1[11]*in2[ 2] + in1[16]*in2[ 3] + in1[21]*in2[ 4];
      tem[ 6]=in1[ 1]*in2[ 5] + in1[ 6]*in2[ 6] + in1[11]*in2[ 7] + in1[16]*in2[ 8] + in1[21]*in2[ 9];
      tem[11]=in1[ 1]*in2[10] + in1[ 6]*in2[11] + in1[11]*in2[12] + in1[16]*in2[13] + in1[21]*in2[14];
      tem[16]=in1[ 1]*in2[15] + in1[ 6]*in2[16] + in1[11]*in2[17] + in1[16]*in2[18] + in1[21]*in2[19];
      tem[21]=in1[ 1]*in2[20] + in1[ 6]*in2[21] + in1[11]*in2[22] + in1[16]*in2[23] + in1[21]*in2[24];
      tem[ 2]=in1[ 2]*in2[ 0] + in1[ 7]*in2[ 1] + in1[12]*in2[ 2] + in1[17]*in2[ 3] + in1[22]*in2[ 4];
      tem[ 7]=in1[ 2]*in2[ 5] + in1[ 7]*in2[ 6] + in1[12]*in2[ 7] + in1[17]*in2[ 8] + in1[22]*in2[ 9];
      tem[12]=in1[ 2]*in2[10] + in1[ 7]*in2[11] + in1[12]*in2[12] + in1[17]*in2[13] + in1[22]*in2[14];
      tem[17]=in1[ 2]*in2[15] + in1[ 7]*in2[16] + in1[12]*in2[17] + in1[17]*in2[18] + in1[22]*in2[19];
      tem[22]=in1[ 2]*in2[20] + in1[ 7]*in2[21] + in1[12]*in2[22] + in1[17]*in2[23] + in1[22]*in2[24];
      tem[ 3]=in1[ 3]*in2[ 0] + in1[ 8]*in2[ 1] + in1[13]*in2[ 2] + in1[18]*in2[ 3] + in1[23]*in2[ 4];
      tem[ 8]=in1[ 3]*in2[ 5] + in1[ 8]*in2[ 6] + in1[13]*in2[ 7] + in1[18]*in2[ 8] + in1[23]*in2[ 9];
      tem[13]=in1[ 3]*in2[10] + in1[ 8]*in2[11] + in1[13]*in2[12] + in1[18]*in2[13] + in1[23]*in2[14];
      tem[18]=in1[ 3]*in2[15] + in1[ 8]*in2[16] + in1[13]*in2[17] + in1[18]*in2[18] + in1[23]*in2[19];
      tem[23]=in1[ 3]*in2[20] + in1[ 8]*in2[21] + in1[13]*in2[22] + in1[18]*in2[23] + in1[23]*in2[24];
      tem[ 4]=in1[ 4]*in2[ 0] + in1[ 9]*in2[ 1] + in1[14]*in2[ 2] + in1[19]*in2[ 3] + in1[24]*in2[ 4];
      tem[ 9]=in1[ 4]*in2[ 5] + in1[ 9]*in2[ 6] + in1[14]*in2[ 7] + in1[19]*in2[ 8] + in1[24]*in2[ 9];
      tem[14]=in1[ 4]*in2[10] + in1[ 9]*in2[11] + in1[14]*in2[12] + in1[19]*in2[13] + in1[24]*in2[14];
      tem[19]=in1[ 4]*in2[15] + in1[ 9]*in2[16] + in1[14]*in2[17] + in1[19]*in2[18] + in1[24]*in2[19];
      tem[24]=in1[ 4]*in2[20] + in1[ 9]*in2[21] + in1[14]*in2[22] + in1[19]*in2[23] + in1[24]*in2[24];
      out[ 0]=tem[ 0];
      out[ 1]=tem[ 1];
      out[ 2]=tem[ 2];
      out[ 3]=tem[ 3];
      out[ 4]=tem[ 4];
      out[ 5]=tem[ 5];
      out[ 6]=tem[ 6];
      out[ 7]=tem[ 7];
      out[ 8]=tem[ 8];
      out[ 9]=tem[ 9];
      out[10]=tem[10];
      out[11]=tem[11];
      out[12]=tem[12];
      out[13]=tem[13];
      out[14]=tem[14];
      out[15]=tem[15];
      out[16]=tem[16];
      out[17]=tem[17];
      out[18]=tem[18];
      out[19]=tem[19];
      out[20]=tem[20];
      out[21]=tem[21];
      out[22]=tem[22];
      out[23]=tem[23];
      out[24]=tem[24];
}

void jam441(double * out,double * in1,double * in2)
{
      double tem[ 4];
      tem[ 0]=in1[ 0]*in2[ 0] + in1[ 4]*in2[ 1] + in1[ 8]*in2[ 2] + in1[12]*in2[ 3];
      tem[ 1]=in1[ 1]*in2[ 0] + in1[ 5]*in2[ 1] + in1[ 9]*in2[ 2] + in1[13]*in2[ 3];
      tem[ 2]=in1[ 2]*in2[ 0] + in1[ 6]*in2[ 1] + in1[10]*in2[ 2] + in1[14]*in2[ 3];
      tem[ 3]=in1[ 3]*in2[ 0] + in1[ 7]*in2[ 1] + in1[11]*in2[ 2] + in1[15]*in2[ 3];
      out[ 0]=tem[ 0];
      out[ 1]=tem[ 1];
      out[ 2]=tem[ 2];
      out[ 3]=tem[ 3];
}
