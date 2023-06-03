/* main.pl */

/* допоміжні функції для обчислень */
abs_FPMIN(X, FPMIN, X1) :- (abs(X) < FPMIN) -> X1 is FPMIN ; X1 is X.

/* підрахунок H */
inc_beta_frac_h(M, MAXIT, H, HRES) :- M > MAXIT -> HRES is 1.0 / 0 ; HRES is H.

/* запуск нового циклу */
launch_new_loop(X, A, B, MAXIT, ESP, FPMIN, M, C, QAB, QAP, QAM, D, H, HRES) :-
    M1 is M + 1,
    M1 < (MAXIT + 1) ->
    inc_beta_loop(X, A, B, MAXIT, ESP, FPMIN, M1, C, QAB, QAP, QAM, D, H, HRES1),
    HRES is HRES1
    ;
    inc_beta_frac_h(M, MAXIT, H, HRES1),
    HRES is HRES1.

/* вихід із циклу */
break_loop(X, A, B, MAXIT, ESP, FPMIN, M, C, QAB, QAP, QAM, D, H, DEL, HRES) :-
    abs(DEL - 1.0) < ESP ->
    inc_beta_frac_h(M, MAXIT, H, HRES1),
    HRES is HRES1
    ;
    launch_new_loop(X, A, B, MAXIT, ESP, FPMIN, M, C, QAB, QAP, QAM, D, H, HRES1),
    HRES is HRES1.

/* тіло циклу, аналогічне функції IncompleteBetaContinuousFraction на мові C */
inc_beta_loop(X, A, B, MAXIT, ESP, FPMIN, M, C, QAB, QAP, QAM, D, H, HRES) :-
    M2 is M * 2,
    AA1 is M * (B - M) * X / ((QAM + M2) * (A + M2)),
    D1 is 1.0 + AA1 * D,
    abs_FPMIN(D1, FPMIN, D_abs1),
    C1 is 1.0 + AA1 / C,
    abs_FPMIN(C1, FPMIN, C_abs1),
    D_rev1 is 1.0 / D_abs1,
    H1 is H * D_rev1 * C_abs1,
    AA2 is (-(A + M) * (QAB + M) * X / ((A + M2) * (QAP + M2))),
    D2 is 1.0 + AA2 * D_rev1,
    abs_FPMIN(D2, FPMIN, D_abs2),
    C2 is 1.0 + AA2 / C_abs1,
    abs_FPMIN(C2, FPMIN, C_abs2),
    D_rev2 is 1.0 / D_abs2,
    DEL is D_rev2 * C_abs2,
    H2 is DEL * H1,
    break_loop(X, A, B, MAXIT, ESP, FPMIN, M, C_abs2, QAB, QAP, QAM, D_rev2, H2, DEL, HRES1),
    HRES is HRES1.

/* передання констант IncompleteBetaContinuousFraction */
inc_beta_frac(X, A, B, H) :-
    QAB is A + B,
    QAP is A + 1.0,
    QAM is A - 1.0,
    D is (1 - QAB * X / QAP),
    abs_FPMIN(D, 1.0e-30, D_abs),
    D1 is 1.0 / D_abs,
    inc_beta_loop(X, A, B, 200, 3.0e-7, 1.0e-30, 1, 1.0, QAB, QAP, QAM, D1, D1, HRES),
    H is HRES.

/* допоміжні вирази для обчислення неповної бета функції */
inc_beta(X, A, B, BETA) :- BETA is e**(lgamma(A + B) - lgamma(A) - lgamma(B) + A * log(X) + B * log(1 - X)).
inc_beta_deter(A, B, DETER) :- DETER is (A + 1) / (A + B + 2).

/* вибір можливого варіанту обчислень неповної бета функції */
inc_beta_branches(X, A, B, BETA, DETER, INCOMPLETE) :-
    X < DETER ->
    inc_beta_frac(X, A, B, H),
    INCOMPLETE is BETA * H / A
    ;
    inc_beta_frac(1 - X, B, A, H),
    INCOMPLETE is 1 - BETA * H / B.

/* підрахунок неповної бета функції */
incomplete_beta_fn(X, A, B, INCOMPLETE) :-
    X < 0 ; X > 1 ->
    INCOMPLETE is 1.0 / 0.0
    ;
    inc_beta(X, A, B, BETA),
    inc_beta_deter(A, B, DETER),
    inc_beta_branches(X, A, B, BETA, DETER, INCOMPLETE1),
    INCOMPLETE is INCOMPLETE1.

/* підрахунок повної бета функції */
complete_beta_fn(A, B, BETA) :- BETA is e**(lgamma(A) + lgamma(B) - lgamma(A + B)).

/* обчислення x* */
count_z(X, MU, SIGMA, Z) :- Z is (X - MU) / SIGMA.
count_xi(X, D1, D2, XI) :- count_z(X, 6.0, 7.0, Z), XI is ((D2 * e**(2*Z)) / (D1 + D2 * e**(2*Z))).

/* F = Ix* */
z_fisher_dist_fn(X, RES) :-
    count_xi(X, 3.0, 2.0, XI),
    incomplete_beta_fn(XI, 3.0 / 2, 2.0 / 2, INCOMPLETE),
    complete_beta_fn(3.0 / 2, 2.0 / 2, BETA),
    RES is INCOMPLETE / BETA.

/* P(a, b) = P(b) - P(a) */
z_fisher_dist_p(A, B, P) :-
    z_fisher_dist_fn(B, F2),
    z_fisher_dist_fn(A, F1),
    P is F2 - F1.

/* допоміжна функція для обчислення ймовірностей */
get_intervals_probability(0, _, _, _, []).
get_intervals_probability(I, Ind, M, Intervals, Probabilities) :-
    I1 is I - 1,
    nth0(Ind, Intervals, Element),
    nth0(0, Element, A),
    nth0(1, Element, B),
    z_fisher_dist_p(A, B, P),
    Ind1 is Ind + 1,
    Probabilities = [P|T],
    get_intervals_probability(I1, Ind1, M, Intervals, T).

/* повернення списку сумарних ймовірностей */
get_probability_list(0, _, _, _, []).
get_probability_list(I, Ind, M, Intervals, PL) :-
    I1 is I - 1,
    Ind1 is Ind + 1,
    nth0(Ind, Intervals, Element),
    get_intervals_probability(M, 0, M, Element, P),
    sumlist(P, Cur_prob),
    PL = [Cur_prob|T],
    get_probability_list(I1, Ind1, M, Intervals, T).

/* допоміжні вирази для роботи зі списками */
copy_list([], []).
copy_list([H|T1], [H|T2]) :- copy_list(T1, T2).

list_size([], 0).
list_size([_|T], Size + 1) :- list_size(T, Size).

replace(I, List, Element, Rep) :-
    nth0(I, List, _, T),
    nth0(I, Rep, Element, T).

index_of([Element|_], Element, 0).
index_of([_|Tail], Element, Index):-
    index_of(Tail, Element, Index1),
    Index is Index1 + 1.

/* допоміжні функція для вибору числа за індексом */
first_interval_border(Arr, Indexes, I, A) :-
    nth0(I, Indexes, Ind),
    nth0(Ind, Arr, Element),
    A is Element.
second_interval_var1(Arr, N, B) :-
    Ind is N - 1,
    nth0(Ind, Arr, Element),
    B is Element.
second_interval_var2(Arr, Indexes, I, B) :-
    nth0(I, Indexes, Ind),
    nth0(Ind, Arr, Element),
    B is Element - 1.
second_interval_border(Arr, Indexes, I, N, M, B) :-
    I < M  ->
    second_interval_var2(Arr, Indexes, I, B1),
    B is B1
    ;
    second_interval_var1(Arr, N, B1),
    B is B1.

/* допоміжна функція для перетворення на інтервали */
translate_to_intervals(0, _, _, _, _, _, []).
translate_to_intervals(I, Ind, Arr, Indexes, N, M, List) :-
    I1 is I - 1,
    Ind1 is Ind + 1,
    first_interval_border(Arr, Indexes, Ind, A),
    second_interval_border(Arr, Indexes, Ind1, N, M, B),
    Interval = [A, B],
    List = [Interval|T],
    translate_to_intervals(I1, Ind1, Arr, Indexes, N, M, T).

/* перетворення інтервалів індексів на інтервали */
indexes_to_intervals(0, _, _, _, _, _, []).
indexes_to_intervals(I, Ind, N, M, Arr, Interval_indexes, Intervals) :-
    I1 is I - 1,
    Ind1 is Ind + 1,
    nth0(Ind, Interval_indexes, Element),
    translate_to_intervals(M, 0, Arr, Element, N, M, List),
    Intervals = [List|T],
    indexes_to_intervals(I1, Ind1, N, M, Arr, Interval_indexes, T).

/* допоміжний вираз для переведення  */
translate_sizes_to_index(0, _, _, _, []).
translate_sizes_to_index(I, Ind, SUM, Sizes_list, Indexes) :-
    I1 is I - 1,
    Indexes = [SUM|T],
    nth0(Ind, Sizes_list, Element),
    SUM1 is SUM + Element,
    Ind1 is Ind + 1,
    translate_sizes_to_index(I1, Ind1, SUM1, Sizes_list, T).

/* повертає переведені в індекси розміри інтервалів */
sizes_to_index_list(0, _, _, _, []).
sizes_to_index_list(I, M, Ind, Interval_sizes, Indexes) :-
    I1 is I - 1,
    add_to_interval_list(M, Ind, Interval_sizes, List),
    translate_sizes_to_index(M, 0, 0, List, Index_list),
    Indexes = [Index_list|T],
    Ind1 is Ind + M,
    sizes_to_index_list(I1, M, Ind1, Interval_sizes, T).

/* допоміжна функція для додавання числа до списку */
add_to_interval_list(0, _, _, []).
add_to_interval_list(I, Ind, Cur_intervals, List) :-
    I1 is I - 1,
    nth0(Ind, Cur_intervals, Element),
    Ind1 is Ind + 1,
    List = [Element|T],
    add_to_interval_list(I1, Ind1, Cur_intervals, T).

/* допоміжна функція для перебору інтервалів */
get_interval_subdiv_sizes(0, _, _, _, _, []).
get_interval_subdiv_sizes(I, Ind, M, Size_index, Cur_intervals, Acc_list) :-
    I1 is I - 1,
    add_to_interval_list(M, Ind, Cur_intervals, List),
    nth0(Size_index, List, Element),
    Iters is Element - 2,
    get_cur_interval_sizes(Iters, Size_index, List, Cur_list),
    Ind1 is Ind + M,
    Acc_list = [Cur_list|T],
    get_interval_subdiv_sizes(I1, Ind1, M, Size_index, Cur_intervals, T).

/* повертає список можливих інтервалів */
get_possible_interval_sizes(0, _, _, _, []).
get_possible_interval_sizes(I, M, Size_index, Cur_intervals, Intervals) :-
    I1 is I - 1,
    flatten(Cur_intervals, Flat_intervals),
    list_size(Flat_intervals, Size),
    Iters is Size / M,
    get_interval_subdiv_sizes(Iters, 0, M, Size_index, Flat_intervals, List),
    Intervals = [List|T],
    Size_index1 is Size_index + 1,
    get_possible_interval_sizes(I1, M, Size_index1, List, T).

/* відповідає тілу цикла GetPossibleIntervalIndexes в C */
get_cur_interval_sizes(-1, _, _, []).
get_cur_interval_sizes(I, S, Intervals, L):-
    S1 is S + 1,
    L = [Intervals|T],
    I1 is I - 1,
    nth0(S, Intervals, Current),
    nth0(S1, Intervals, Next),
    Current1 is Current - 1,
    Next1 is Next + 1,
    replace(S, Intervals, Current1, Intervals1),
    replace(S1, Intervals1, Next1, Intervals2),
    get_cur_interval_sizes(I1, S, Intervals2, T).

/*відбір початкових можливих інтервалів*/
get_cur_interval_sizes_iters(S, MIN_INTERVAL_NUMBER, Intervals, L) :-
    nth0(S, Intervals, First),
    I is First - MIN_INTERVAL_NUMBER,
    get_cur_interval_sizes(I, S, Intervals, L).

get_init_possible_interval_sizes(Init_intervals, MIN_INTERVAL_NUMBER, Current_interval_sizes) :-
    copy_list(Init_intervals, Cur_intervals),
    get_cur_interval_sizes_iters(0, MIN_INTERVAL_NUMBER, Cur_intervals, Current_interval_sizes).

/* початкові розмірності інтервалів */
count_init_interval_size(I, N, M, MIN_INTERVAL_SIZE, RES) :-
    I == 0 -> RES is N - (M - 1) * MIN_INTERVAL_SIZE ; RES is MIN_INTERVAL_SIZE.
get_init_interval_sizes(0, _, _, _, []).
get_init_interval_sizes(I, N, M, MIN_INTERVAL_SIZE, L) :-
    I > 0,
    I1 is I - 1,
    count_init_interval_size(I1, N, M, MIN_INTERVAL_SIZE, RES),
    L = [RES|T],
    get_init_interval_sizes(I1, N, M, MIN_INTERVAL_SIZE, T).

/* допоміжна функція для перекладу числа на літеру */
translate_number_to_letter(0, _, _, _, []).
translate_number_to_letter(I, Alphabet, Num, Intervals, Letter) :-
    I1 is I - 1,
    nth0(I1, Intervals, Interval),
    nth0(0, Interval, A),
    A =< Num ->
    I_next is 0,
    nth0(I1, Alphabet, Element),
    Letter = [Element|T],
    translate_number_to_letter(I_next, Alphabet, Num, Intervals, T)
    ;
    I_next is I - 1,
    Letter = [[]|T],
    translate_number_to_letter(I_next, Alphabet, Num, Intervals, T).

/* переклад чисел на лінгвістичний ряд */
translate_to_ling_chain(0, _, _, _, _, _, []).
translate_to_ling_chain(I, Ind, M, Alphabet, Nums, Intervals, Ling) :-
    I1 is I - 1,
    Ind1 is Ind + 1,
    nth0(Ind, Nums, Num),
    translate_number_to_letter(M, Alphabet, Num, Intervals, Letter),
    flatten(Letter, L),
    Ling = [L|T],
    translate_to_ling_chain(I1, Ind1, M, Alphabet, Nums, Intervals, T).

lab3_main() :-
    Numbers = [13, 1, 2, 15, -10, 20, 3, 12, 4, 5, 6, 26, 7, 8],
    Alphabet = ["A", "B", "C", "D"],
    list_size(Numbers, N),
    list_size(Alphabet, M),

    statistics(walltime, _),

    sort(Numbers, Sorted),

    get_init_interval_sizes(M, N, M, 2, R_Init_interval_sizes),
    reverse(R_Init_interval_sizes, Init_interval_sizes),
    get_init_possible_interval_sizes(Init_interval_sizes, 2, Init_possible_interval_sizes),

    Iters is M - 2,
    get_possible_interval_sizes(Iters, M, 1, Init_possible_interval_sizes, Possible_interval_sizes),
    flatten(Possible_interval_sizes, Flat_interval_sizes),
    list_size(Flat_interval_sizes, S),
    Ind_iters is S / M,
    sizes_to_index_list(Ind_iters, M, 0, Flat_interval_sizes, Interval_indexes),
    list_size(Interval_indexes, IS),
    indexes_to_intervals(IS, 0, N, M, Sorted, Interval_indexes, Intervals),

    list_size(Intervals, II),
    get_probability_list(II, 0, M, Intervals, Probabilities),

    max_list(Probabilities, Max_e),
    index_of(Probabilities, Max_e, Max_ind),
    nth0(Max_ind, Intervals, Probable_intervals),

    MN is M + 0,
    translate_to_ling_chain(N, 0, MN, Alphabet, Numbers, Probable_intervals, Ling),
    flatten(Ling, Chain),
    writeln("Linguistic chain"),
    writeln(Chain),

    statistics(walltime, [_ | Time]),
    write("Execution time "), write(Time), writeln(" ms").
