/* підрахунок початкових станів змінних IncompleteBetaContinuousFraction */
incBetaQAB(A, B, QAB) :- QAB is A + B.
incBetaQAP(A, QAP) :- QAP is A + 1.0.
incBetaQAM(A, QAM) :- QAM is A - 1.0.
incBetaD(X, QAB, QAP, D) :- D is (1 - QAB * X / QAP).

/* допоміжні функції для обчислень */
absFPMIN(X, FPMIN, X1) :- (abs(X) < FPMIN) -> X1 is FPMIN ; X1 is X.
reverseX(X, X1) :- X1 is 1 / X.
mulM2(M, M2) :- M2 is M * 2.
countAA1(X, A, B, M, M2, QAM, AA) :- AA is M * (B - M)*X / ((QAM + M2)*(A + M2)).
countD1(D, AA, D1) :- D1 is 1.0 + AA * D.
countC1(C, AA, C1) :- C1 is 1.0 + AA / C.
countH1(H, D, C, H1) :- H1 is H * D * C.
countDEL(D, C, DEL) :- DEL is D * C.
countH2(H, DEL, H2) :- H2 is DEL * H.
countAA2(X, A, M, M2, QAB, QAP, AA2) :- AA2 is (-(A + M)*(QAB + M) * X / ((A + M2) * (QAP + M2))).
increment(X, X1) :- X1 is X + 1.

/* підрахунок H */
incBetaFracH(M, MAXIT, H, HRES) :-
    M > MAXIT ->
    HRES is 1.0 / 0
    ;
    HRES is H.

/* запуск нового циклу */
launchNewLoop(X, A, B, MAXIT, ESP, FPMIN, M, C, QAB, QAP, QAM, D, H, HRES) :-
    increment(M, M1),
    M1 < (MAXIT + 1) ->
    incBetaLoop(X, A, B, MAXIT, ESP, FPMIN, M1, C, QAB, QAP, QAM, D, H, HRES1),
    HRES is HRES1
    ;
    incBetaFracH(M, MAXIT, H, HRES1),
    HRES is HRES1.

/* вихід із циклу */
breakLoop(X, A, B, MAXIT, ESP, FPMIN, M, C, QAB, QAP, QAM, D, H, DEL, HRES) :-
    abs(DEL - 1.0) < ESP ->
    incBetaFracH(M, MAXIT, H, HRES1),
    HRES is HRES1
    ;
    launchNewLoop(X, A, B, MAXIT, ESP, FPMIN, M, C, QAB, QAP, QAM, D, H, HRES1),
    HRES is HRES1.

/* тіло циклу, аналогічне функції IncompleteBetaContinuousFraction на мові C */
incBetaLoop(X, A, B, MAXIT, ESP, FPMIN, M, C, QAB, QAP, QAM, D, H, HRES) :-
    mulM2(M, M2),
    countAA1(X, A, B, M, M2, QAM, AA1),

    countD1(D, AA1, D1),
    absFPMIN(D1, FPMIN, Dabs1),

    countC1(C, AA1, C1),
    absFPMIN(C1, FPMIN, Cabs1),

    reverseX(Dabs1, Drev1),

    countH1(H, Drev1, Cabs1, H1),
    countAA2(X, A, M, M2, QAB, QAP, AA2),

    countD1(Drev1, AA2, D2),
    absFPMIN(D2, FPMIN, Dabs2),

    countC1(Cabs1, AA2, C2),
    absFPMIN(C2, FPMIN, Cabs2),

    reverseX(Dabs2, Drev2),

    countDEL(Drev2, Cabs2, DEL),
    countH2(H1, DEL, H2),

    breakLoop(X, A, B, MAXIT, ESP, FPMIN, M, Cabs2, QAB, QAP, QAM, Drev2, H2, DEL, HRES1),
    HRES is HRES1.

/* передання констант IncompleteBetaContinuousFraction */
incBetaFrac(X, A, B, H) :-
    incBetaQAB(A, B, QAB),
    incBetaQAP(A, QAP),
    incBetaQAM(A, QAM),
    incBetaD(X, QAB, QAP, D),
    absFPMIN(D, 1.0e-30, Dabs),
    reverseX(Dabs, D1),
    incBetaLoop(X, A, B, 200, 3.0e-7, 1.0e-30, 1, 1.0, QAB, QAP, QAM, D1, D1, HRES),
    H is HRES.

/* допоміжні вирази для обчислення неповної бета функції */
incBeta(X, A, B, BETA) :- BETA is e**(lgamma(A + B) - lgamma(A) - lgamma(B) + A * log(X) + B * log(1 - X)).
incBetaDeter(A, B, DETER) :- DETER is (A + 1) / (A + B + 2).


/* вибір можливого варіанту обчислень неповної бета функції */
incBetaBranches(X, A, B, BETA, DETER, INCOMPLETE) :-
    X < DETER ->
    incBetaFrac(X, A, B, H),
    INCOMPLETE is BETA * H / A
    ;
    incBetaFrac(1 - X, B, A, H),
    INCOMPLETE is 1 - BETA * H / B.

/* підрахунок неповної бета функції */
incompleteBetaFn(X, A, B, INCOMPLETE) :-
    X < 0 ; X > 1 ->
    INCOMPLETE is 1.0 / 0.0
    ;
    incBeta(X, A, B, BETA),
    incBetaDeter(A, B, DETER),
    incBetaBranches(X, A, B, BETA, DETER, INCOMPLETE1),
    INCOMPLETE is INCOMPLETE1.

/* підрахунок повної бета функції */
completeBetaFn(A, B, BETA) :- BETA is e**(lgamma(A) + lgamma(B) - lgamma(A + B)).

/* обчислення x* */
countZ(X, MU, SIGMA, Z) :- Z is (X - MU) / SIGMA.
countXI(X, D1, D2, XI) :- countZ(X, 6.0, 7.0, Z), XI is ((D2 * e**(2*Z)) / (D1 + D2 * e**(2*Z))).

/* F = Ix* */
zFisherDistFn(X, RES) :-
    countXI(X, 3.0, 2.0, XI),
    incompleteBetaFn(XI, 3.0 / 2, 2.0 / 2, INCOMPLETE),
    completeBetaFn(3.0 / 2, 2.0 / 2, BETA),
    RES is INCOMPLETE / BETA.

/* P(a, b) = P(b) - P(a) */
zFisherDistP(A, B, P) :-
    zFisherDistFn(B, F2),
    zFisherDistFn(A, F1),
    P is F2 - F1.

