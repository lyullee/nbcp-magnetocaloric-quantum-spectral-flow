# XTRG 프로토타입 결과 사용 구분

## 결론에 사용할 기준 파일

- `xtrg_quimb_n9_purified_stabilized_convergence.json`
  - 올바른 MPO 인덱스 방향을 사용한다.
  - `rho=X†X` 양성 표현을 사용한다.
  - D=16, 24, 32의 같은 조건 비교와 `STOP_FOR_BULK` 판정을 포함한다.
- `runtime_tensor_network_manifest_20261011.json`
  - Octave 파일 해시, Windows 차단 원인, Python 런타임을 기록한다.
- 그림: `../figures/xtrg_quimb_n9_feasibility.png`

## 진단용으로만 보존한 파일

- 이름에 `INVALID_INDEX_ORIENTATION`가 있는 파일은 `.H`의 인덱스 방향을 잘못 해석한 실행이다. 과학 결과로 사용하지 않는다.
- 이름에 `UNSTABILIZED`가 있는 파일은 정화 인자에서 압축 비에르미트 오차를 매 단계 제거하지 않은 실행이다.
- `xtrg_quimb_n9_validation_D32_raw_square.json`과 `xtrg_quimb_n9_validation.json`은 직접 밀도 MPO 압축의 실패 진단이다.

이전 파일을 삭제하지 않은 이유는 오류 발견 과정과 중단결정을 감사 가능하게 남기기 위해서다.
