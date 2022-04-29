/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 22-03-2021
Language/ver     : C++ in MSVS2019

Description      : NM_MainTemplate.cpp
-------------------------------------------------------------------------------*/

#define Assignment   4      // enter your assignment number
#define eval      0         // set

#include "myNM.h"

int main(int argc, char* argv[])
{
    /*    [¡Ø DO NOT EDIT IT !!!]   Resources file path setting for evaluation   */
    std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";

#if eval
    path += "eval/";
#endif

    /*==========================================================================*/
    /*               Variables declaration & initialization               */
    /*--------------------------------------------------------------------------*/
    /*   - You can change the variable names                            */
    /*   - However, you must use the specified file name                  */
    /*      : For each assignment, the file name will be notified on HISNET      */
    /*==========================================================================*/

    Matrix matA = txt2Mat(path, "matA");
    Matrix matC = txt2Mat(path, "matC");

    ///*==========================================================================*/
    ///*               Apply your numerical method algorithm                 */
    ///*==========================================================================*/

    // Condition Number Matrix
    Matrix transposeA = createMat(matA.cols, matA.rows);
    Matrix ATA = createMat(matA.cols, matA.cols);
    Matrix tempA = createMat(matA.cols, matA.cols);
    Matrix eigenATA = createMat(ATA.rows, ATA.cols);

    Matrix transposeC = createMat(matC.cols, matC.rows);
    Matrix CTC = createMat(matC.cols, matC.cols);
    Matrix tempC = createMat(matC.cols, matC.cols);
    Matrix eigenCTC = createMat(CTC.rows, CTC.cols);

    double Cond_NumberA = 0;
    double Cond_NumberC = 0;

    /*==========================================================================*/
    /*                             Input Matrix A                              */
    /*==========================================================================*/
    printf("\n-------------------------------------------------------------");
    printf("\n                      Input Matrix A                         ");
    printf("\n-------------------------------------------------------------\n");
    printMat(matA, "[ matA ]");


    /*==========================================================================*/
    /*                            Input Matrix C                              */
    /*==========================================================================*/
    printf("\n-------------------------------------------------------------");
    printf("\n                      Input Matrix C                         ");
    printf("\n-------------------------------------------------------------\n");
    printMat(matC, "[ matC ]");


    /*==========================================================================*/
    /*                            Find ATA                                      */
    /*==========================================================================*/
    printf("\n-------------------------------------------------------------");
    printf("\n                   Matrix     A T A                          ");
    printf("\n-------------------------------------------------------------\n");
    transposeA = transpose(matA);
    multiply_Mat(transposeA, matA, tempA);
    ATA = copyMat(tempA);
    printMat(ATA, "[ ATA ]");


    /*==========================================================================*/
    /*                            Find CTC                                      */
    /*==========================================================================*/
    printf("\n-------------------------------------------------------------");
    printf("\n                   Matrix     C T C                          ");
    printf("\n-------------------------------------------------------------\n");
    transposeC = transpose(matC);
    multiply_Mat(transposeC, matC, tempC);
    CTC = copyMat(tempC);
    printMat(CTC, "[ CTC ]");


    /*==========================================================================*/
    /*                            Print EigenValue                          */
    /*==========================================================================*/
    // QR Decomposition Matrix Setting
    Matrix ATAmatH = createMat(ATA.rows, ATA.cols);
    Matrix ATAmatQ = createMat(ATA.rows, ATA.cols);
    Matrix ATAmatR = createMat(ATA.rows, ATA.cols);
    printf("\n--------------------------------------------------------------");
    printf("\n                EigenValues  of  A T A                        ");
    printf("\n--------------------------------------------------------------\n");
    eigenATA = eig(ATA, ATAmatH, ATAmatQ, ATAmatR);
    printMat(eigenATA, "eigenValue of ATA");


    /*==========================================================================*/
    /*                            Print EigenValue                          */
    /*==========================================================================*/
    // QR Decomposition Matrix Setting
    Matrix CTCmatH = createMat(CTC.rows, CTC.cols);
    Matrix CTCmatQ = createMat(CTC.rows, CTC.cols);
    Matrix CTCmatR = createMat(CTC.rows, CTC.cols);
    printf("\n--------------------------------------------------------------");
    printf("\n                EigenValues  of  C T C                        ");
    printf("\n--------------------------------------------------------------\n");
    eigenCTC = eig(CTC, CTCmatH, CTCmatQ, CTCmatR);
    printMat(eigenCTC, "eigenValue of CTC");


    /*==========================================================================*/
    /*                  Find Maximum and Minimum EigenValues                    */
    /*==========================================================================*/
    double maxeig_ATA = 0;
    double maxeig_CTC = 0;
    double mineig_ATA = 0;
    double mineig_CTC = 0;
    maxeig_ATA = max(eigenATA);
    maxeig_CTC = max(eigenCTC);
    mineig_ATA = min(eigenATA);
    mineig_CTC = min(eigenCTC);

    printf("\n--------------------------------------------------------------");
    printf("\n                Max | Min Eigenvalue of A T A                 ");
    printf("\n--------------------------------------------------------------\n");
    printf("Max = %f\t", maxeig_ATA);
    printf("Min = %f\t", mineig_ATA);
    printf("\n\n");

    printf("\n--------------------------------------------------------------");
    printf("\n                Max | Min Eigenvalue of C T C                 ");
    printf("\n--------------------------------------------------------------\n");
    printf("Max = %f\t", maxeig_CTC);
    printf("Min = %f\t", mineig_CTC);
    printf("\n\n");


    /*==========================================================================*/
    /*                            Condition Number                              */
    /*==========================================================================*/
    printf("\n-------------------------------------------------------------");
    printf("\n                 Condition Number of Mat A                   ");
    printf("\n-------------------------------------------------------------\n");
    Cond_NumberA = condNumber(matA);
    printf("%f\t", Cond_NumberA);
    printf("\n\n");

    printf("\n-------------------------------------------------------------");
    printf("\n                 Condition Number of Mat C                   ");
    printf("\n-------------------------------------------------------------\n");
    Cond_NumberC = condNumber(matC);
    printf("%f\t", Cond_NumberC);
    printf("\n\n");


    /*==========================================================================*/
    /*                       Deallocate memory                      */
    /*==========================================================================*/
    freeMat(matA);          freeMat(matC);          
    freeMat(transposeA);    freeMat(ATA);           freeMat(tempA);         freeMat(eigenATA);
    freeMat(transposeC);    freeMat(CTC);           freeMat(tempC);         freeMat(eigenCTC);
    freeMat(ATAmatH);       freeMat(ATAmatQ);       freeMat(ATAmatR);
    freeMat(CTCmatH);       freeMat(CTCmatQ);       freeMat(CTCmatR);

    system("pause");
    return 0;
}


