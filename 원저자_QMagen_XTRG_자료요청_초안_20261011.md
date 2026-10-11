# QMagen/XTRG 실행환경 및 수치자료 요청 초안

보내기 전 수신자와 문구를 사용자가 검토해야 한다. 현재 공개 QMagen README에 적힌 연락처는 Yuan Gao의 `17231064@buaa.edu.cn`이다.

## 영문 제목

Request for numerical XTRG data and reproducibility details for Na2BaCo(PO4)2

## 영문 본문

Dear Dr. Gao,

I am studying how weak realistic interactions, including bond-dependent exchange and long-range dipolar coupling, affect the sub-kelvin magnetocaloric operating window of Na2BaCo(PO4)2. I have reproduced the bulk M0 operating-field window by digitizing Fig. 4 of your 2022 npj Quantum Materials paper, and I found the public QMagen repository and its MATLAB XTRG implementation.

To make the comparison quantitative and reproducible, could you please share or clarify the following items?

1. The numerical arrays used for Fig. 4a and Fig. 4b: temperature, magnetic field, C_m/T, and entropy.
2. The exact QMagen commit, MATLAB release, operating system, and required toolboxes used for the YC6×9 XTRG calculation.
3. The initial tau or beta, SETTN order, bond-dimension schedule, truncation cutoff, number of XTRG iterations, and convergence criteria.
4. The site ordering and boundary convention used for the YC6×9 triangular cylinder.
5. Whether the public interaction-map/MPO builder has been tested with complex bond-dependent SxSy/SySx terms or long-range pair interactions.
6. The intended software license for the public QMagen repository, which currently does not appear to contain a LICENSE file.

I will cite the original paper and repository, preserve all provenance, and use any shared numerical data only under the conditions you specify. A machine-readable table for Fig. 4 would already be very helpful even if the original run environment is no longer available.

Sincerely,

Ugwiyeon Lee  
Korea Gas Safety Corporation  
ORCID: 0009-0008-8976-5363  
lyullee@kgs.or.kr

## 요청자료 체크리스트

- `Fig4a.csv`: `T_K, B_T, Cm_over_T_J_mol_K2`
- `Fig4b.csv`: `T_K, B_T, entropy_J_mol_K`
- 실행환경: OS, MATLAB 버전, toolbox, CPU/RAM
- 실행입력: lattice/site ordering, BC, Hamiltonian parameters, field convention
- 수렴입력: `tau`, SETTN order, `D` schedule, cutoff, iteration count
- 원시 출력: `beta, lnZ, U, C, S`, truncation error, wall time
- 저장소 라이선스 또는 재사용 허용범위
