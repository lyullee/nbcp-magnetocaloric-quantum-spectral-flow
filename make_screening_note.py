"""Generate a plain-language record from completed numerical outputs."""
from pathlib import Path
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'results'
summary = json.loads((OUT/'screening_summary.json').read_text(encoding='utf-8'))
weak = json.loads((OUT/'weak_term_sensitivity.json').read_text(encoding='utf-8'))
manifest = json.loads((ROOT/'source_manifest.json').read_text(encoding='utf-8'))
for item in manifest['files']:
    data=(ROOT/'sources'/item['file']).read_bytes()
    assert len(data)==item['bytes']
    assert hashlib.sha256(data).hexdigest()==item['sha256']

audit = summary['source_data_audit']
candidates = audit['measured_field_candidates_100mK']
target = audit['conditional_decision_accuracy_target']
with (OUT/'finite_cluster_cycle_profiles.csv').open(encoding='utf-8-sig',newline='') as f:
    profiles=list(csv.DictReader(f))
matched=[]
for n in [9,12]:
    for model in ['XXZ','Ising']:
        rows=[r for r in profiles if int(r['N'])==n and r['model']==model
              and abs(float(r['T_cold_K'])-.1)<1e-9
              and any(abs(float(r['B_end_T'])-b)<1e-9 for b in [0.,1.3,1.5,1.66])]
        assert len(rows)==4
        best=max(rows,key=lambda r:float(r['Q_entropy_budget_mJ_per_cm3']))
        matched.append({'N':n,'model':model,'best_of_four_candidates_B_T':float(best['B_end_T']),
            'rows':[{k:(float(v) if k not in ['model'] else v) for k,v in r.items()} for r in rows]})
(OUT/'matched_field_decisions.json').write_text(json.dumps(matched,indent=2),encoding='utf-8')

lines=[
'자기냉동 후보 검증: 공개 데이터와 로컬 계산으로 확인한 내용',
'작성일: 2026-10-11',
'상태: 연구문제의 가능성을 검토한 단계. 논문 결과 또는 양자우위 증명은 아님.',
'',
'1. 우리가 해결하려는 공학 문제',
'냉동기는 차가운 곳의 열을 가져와 밖으로 버리는 장치다. 자기냉동기는 자석의 세기를 바꾸어 자성체가 열을 받아들이거나 내놓게 한다.',
'여기서 목표는 100 mK, 즉 0.1 K에서 사용할 냉매의 운전 자기장과 흡수 가능한 열량을 예측하는 것이다. 수소 액화 온도 약 20 K와는 다른 초저온 장치다.',
'냉매 내부의 스핀은 원자에 붙은 작은 양자 자석으로 생각할 수 있다. 이 양자 자석들이 서로 방해하거나 도와주는 방식이 냉각성능에 영향을 줄 수 있다.',
'현재 연구 질문: 실제 재료에서 작은 상호작용까지 포함했을 때, 단순한 모델이 고른 운전조건은 유지되는가? 더 정확한 양자 열상태 계산은 어떤 조건에서 설계에 필요한가?',
'열상태는 정해진 온도에서 여러 에너지 상태가 알맞은 비율로 섞인 상태다. 양자컴퓨터가 맡을 후보 작업은 이 상태의 에너지와 엔트로피를 계산하는 것이다.',
'',
'2. 실제로 완료한 작업',
'Nature 2024 논문이 공개한 원본 데이터 파일 3개, 총 18개 시트를 내려받았다. 파일 크기와 SHA-256 지문을 저장하고 원본이 바뀌지 않았음을 다시 확인했다.',
'원본 실험실 로그 전체를 확보한 것은 아니다. 논문 그림에 사용된 공개 수치 데이터다. 실험값, 보정값, 이론값을 구분하여 사용했다.',
'고전 PC에서 9개와 12개 스핀의 정확한 에너지 상태와 열역학 값을 계산했다. 작은 계의 정확한 답이며 큰 재료의 답으로 확장한 결과는 아니다.',
'양자회로 또는 실제 양자컴퓨터는 실행하지 않았다. 유료 클라우드 계산을 신청하지 않았다.',
f'기본 모델 계산 시간: {summary["elapsed_seconds"]:.2f}초. 약한 상호작용 민감도 계산 시간: {weak["elapsed_seconds"]:.2f}초.',
'검산: 에너지 행렬의 대칭성, 상태 수, 높은 온도의 엔트로피, 큰 자기장의 정렬 상태, 에너지/엔트로피의 온도 미분을 확인했다. 추가 모델은 별도로 만든 두 계산 방식의 고유값을 대조했다.',
'',
'3. 공개 데이터에서 계산한 운전조건 비교',
'공통 출발 상태: 2 K, 4 T. 목표 온도: 0.1 K. 아래는 공개된 네 자기장 후보의 비교다.',
'엔트로피는 재료 내부에 남아 있는 여러 가능한 배열의 정도를 나타내는 값이다. 냉매가 열을 받을 때 그 배열이 더 다양해질 수 있다.',
'계산식: Q = T_cold × (S_cold - S_hot). Q는 냉매 부피당 열량이며, 출발 상태에서 냉각 후 등온 흡열 구간으로 가는 가역 경로가 가능하다는 조건이 붙는다.',
'자기장(T) | 0.1 K의 엔트로피(mJ/cm3/K) | 흡열 구간의 열량 예산(mJ/cm3)',
]
for r in candidates:
    lines.append(f'{r["B_end_T"]:.2f} | {r["S_cold_mJ_per_cm3_K"]:.6f} | {r["Q_entropy_budget_mJ_per_cm3"]:.6f}')
lines += [
'네 후보의 중앙값 비교에서는 1.5 T가 가장 유리하다. 자기장을 완전히 끄는 것이 항상 유리하지 않다는 실제 설계 질문이 된다.',
'이 값은 새로 측정한 냉각성능이 아니다. 공개 엔트로피 자료를 보간하여 계산한 결과다. 냉각 속도, W 단위의 냉각능력, 전기 효율 또는 전체 사이클의 최적조건을 뜻하지 않는다.',
'출처: extended_figure_6.xlsx의 c 시트(A/B: 0 T, C/D: 4 T), extended_figure_4.xlsx의 h/i/j 시트(A: 온도, D: 실험 기반 자기 엔트로피). 보간 범위 안의 값만 사용했다.',
'',
'4. 작은 모델이 실제 설계 선택을 맞히는가?',
]
for r in matched:
    lines.append(f'N={r["N"]}, {r["model"]}: 동일한 네 자기장 중 {r["best_of_four_candidates_B_T"]:.2f} T를 선택.')
lines += [
'XXZ는 스핀의 양자적 교환을 포함한 모델이다. Ising은 그 일부를 지운 비교 모델이다. Ising과의 차이는 물리를 지운 효과이지, 고전 알고리즘과 양자 알고리즘의 계산 성능 차이가 아니다.',
'9개와 12개 스핀은 공개 데이터가 고른 1.5 T를 그대로 재현하지 못했다. 적은 큐비트의 실행 결과를 곧바로 실제 냉동기 설계 결과라고 주장할 근거가 없다는 뜻이다.',
'다음 계산에는 계의 크기, 경계 모양, 상호작용 범위를 늘렸을 때 답이 안정되는지 확인하는 절차가 필요하다.',
'',
'5. 약한 상호작용을 추가한 민감도 계산',
'2026년 공개 원고의 별도 매개변수 세트를 사용했다. 2024년 모델에 서로 다른 매개변수를 섞지 않았다.',
'비교는 (a) 기본 XXZ, (b) 추가적인 결합 방향별 상호작용 Jpp, (c) 가까운 이웃의 쌍극자 상호작용과 추가 Jpp를 함께 넣은 모델이다.',
'쌍극자 상호작용은 작은 자석끼리 주고받는 힘이다. 결합 방향별 상호작용은 두 스핀이 연결된 방향에 따라 달라지는 힘이다.',
'이는 2차원 소형 계의 민감도 검사다. 실제 원고의 멀리 떨어진 스핀 간 상호작용과 층 사이 상호작용을 모두 넣은 3차원 모델은 아니다.',
'아래 변화량은 같은 크기의 기본 XXZ 결과에서 얼마나 달라졌는지 뜻한다.',
'N | 자기장(T) | 추가 모델 | 엔트로피 변화(kB/스핀) | 열량 예산 변화(mJ/cm3)',
]
for c in weak['comparisons']:
    label='추가 Jpp' if c['model']=='XXZ_plus_additional_Jpp' else '가까운 쌍극자+추가 Jpp'
    lines.append(f'{c["N"]} | {c["B_end_T"]:.2f} | {label} | {c["delta_S_cold_kB_per_spin"]:+.6f} | {c["delta_Q_entropy_budget_mJ_per_cm3"]:+.6f}')
lines += [
'추가 Jpp만 넣었을 때 냉각 끝점 엔트로피의 변화는 이번 계산에서 스핀당 약 0.00038 kB 이하로 작았다. 이 항 하나를 넣는 것만으로 강한 공학적 결과를 확보했다고 말할 수 없다.',
'가까운 쌍극자 힘까지 포함한 경우 1.5 T의 열량 예산은 9개 스핀에서 약 +0.204 mJ/cm3, 12개 스핀에서 약 -0.081 mJ/cm3 바뀌었다. 변화 방향부터 일치하지 않는다. 작은 모델에서의 민감도와 유한 크기의 영향이 섞여 있다는 신호다.',
'이 작은 계의 변화량만으로 실제 냉매의 성능 변화를 확정할 수 없다. 크기를 바꾸었을 때 결과가 일치하는지, 긴 거리와 층 사이의 상호작용이 어떤 영향을 주는지 추가 검증해야 한다.',
'',
'6. 계산 정확도를 공학적으로 정하는 방법',
f'1.5 T와 다음 후보 1.66 T의 열량 예산 차이는 {target["Q_gap_mJ_per_cm3"]:.6f} mJ/cm3다.',
f'다른 오차가 없다고 가정하면, 두 후보 각각의 엔트로피 절대오차가 {target["equal_absolute_entropy_error_bound_each_endpoint_mJ_per_cm3_K"]:.6f} mJ/cm3/K보다 작을 때 이 두 후보의 중앙값 순서를 보존할 수 있다.',
f'재료의 자기 이온 밀도로 환산하면 스핀당 {target["equal_absolute_entropy_error_bound_each_endpoint_kB_per_spin"]:.6f} kB다.',
'이는 실험 신뢰구간이 아니다. 공개 자료의 측정/보정 오차, 보간 오차, 모델 매개변수 불확실성이 추가된다. 최적 자기장 전체를 찾아낸 기준도 아니다.',
'따라서 양자 알고리즘의 정확도를 임의로 정할 필요가 없다. 실제 운전조건의 선택이 바뀌지 않도록 필요한 정확도와 측정 횟수를 정할 수 있다.',
'',
'7. 문헌에서 확인한 기회와 경쟁 상대',
'2024년 Nature 연구는 이미 고전 텐서네트워크 계산과 실험으로 이 재료의 큰 자기냉각 효과를 다뤘다. 같은 단순 모델의 냉각 효과를 다시 계산하는 것만으로는 새로운 연구가 되기 어렵다.',
'2026년 공개 원고는 약 1.35-1.65 T 구간에서 단순 XXZ 모델이 설명하지 못하는 새로운 자성 배열을 보고했다. 이 구간은 위의 1.5 T 운전 후보를 포함한다. 공개 원고이므로 심사 완료 논문과 구분해야 한다.',
'그 연구는 약한 쌍극자 힘과 결합 방향별 힘, 층 사이 힘을 넣어 반고전적으로 자성 배열을 설명한다. 반고전은 양자효과의 일부를 근사해서 계산하는 방법이다. 그 배열을 발견했다는 이유만으로 양자컴퓨터가 필요하다고 주장할 수 없다.',
'현재 검토한 범위에서는 그 현실적인 상호작용이 유한 온도의 엔트로피와 냉각 운전조건에 미치는 영향을 양자 열상태 알고리즘과 연결한 완결된 연구는 확인하지 못했다. 신규성은 추가 확인 대상이다.',
'층 사이 힘의 값은 한 가지로 정해지지 않았다. 여러 가능한 값에서도 결론이 유지되는지 검증해야 한다.',
'다른 삼각격자 자성체 K2Co(SeO3)2에서는 계산을 쉽게 만드는 저에너지 모델과 고전 몬테카를로가 실험을 잘 설명한다. 복잡한 자성체라는 이유만으로 고전 계산이 불가능하다고 해서는 안 된다.',
'현재 판단: 공개 데이터 기반의 검증은 시작할 수 있고 실제로 시작했다. 의미 있는 공학적 차이를 얻을 가능성은 아직 검증 중이다. 양자컴퓨터가 고전 최선 방법보다 유리하다는 근거는 아직 확보하지 않았다. 이 재료를 논문 주제로 최종 확정한 단계는 아니다.',
'',
'8. 다음 단계와 중단/전환 기준',
'A. 공개 열역학 자료의 보정/오차를 반영해 운전 자기장 선택의 확실성을 계산한다. 같은 데이터를 맞추는 데 쓴 부분과 예측 검증용 부분을 구분한다.',
'B. 단순 2D 모델 → 결합 방향별 힘 → 긴 거리 쌍극자 힘 → 층 사이 힘을 순서대로 추가한다. 크기와 경계 조건을 바꾸고, 동일한 모델에 대해 텐서네트워크/저에너지 근사/반고전 계산과 비교한다.',
'C. 모델 변화와 양자효과가 설계 선택을 바꾸고, 그 변화가 데이터와 매개변수의 불확실성보다 큰 구간을 찾는다. 이 조건이 충족되지 않으면 이 재료를 주제로 고정하지 않는다.',
'D. 그 구간에서 필요한 열상태 준비, 에너지 측정, 엔트로피 복원 비용을 계산한다. 논리 큐비트, 보조 큐비트, 회로 깊이, 측정 횟수, 온도별 반복 횟수까지 포함한다.',
'E. 작은 계에서 회로 시뮬레이터가 정확한 답과 필요한 설계 정확도를 재현하는지 먼저 확인한다. 이후 무료 QPU에서 결과 해석에 필요한 소수 회로만 실행한다. 작은 QPU 결과를 큰 재료의 성능 실증으로 부르지 않는다.',
'F. 현재 Q 수식의 조건부 열량 예산에서 실제 냉동 사이클로 확장하려면 경로 도달 가능성, 열스위치, 외부 열 유입, 자성체/격자 열교환을 추가해야 한다.',
'',
'근거 문헌과 데이터',
'Nature 2024: https://www.nature.com/articles/s41586-023-06885-w',
'2024 논문 저자 원고/모델: https://arxiv.org/html/2504.11298v1',
'2026 새로운 자성 배열과 확장 모델(공개 원고): https://arxiv.org/html/2607.23173v1',
'고전 저에너지 모델의 성공 사례: https://www.nature.com/articles/s41535-025-00791-2',
'유한 온도의 반고전/양자 보정 경쟁 계산: https://arxiv.org/html/2505.07229v1',
'유효 모델의 적용 범위에 관한 2026 연구: https://doi.org/10.1103/f5xy-9q3m',
'2026 원고가 인용한 ISIS 데이터 DOI: https://doi.org/10.5286/ISIS.E.RB2010721',
'주의: 이 DOI의 카탈로그 제목은 CoTiO3/FeTiO3 실험이다. 해당 자료 안에 NBCP의 필요한 런이 있는지 확인하지 못했으므로 이번 계산에는 쓰지 않았다.',
'원본 파일의 직접 주소와 SHA-256: source_manifest.json',
'재실행 코드: fetch_sources.py, screen_model.py, screen_weak_terms.py, make_screening_note.py',
'계산 결과: results 폴더의 CSV와 JSON',
]
(ROOT/'자기냉동_검증결과_쉽게설명.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8-sig')
print(json.dumps({'report':str(ROOT/'자기냉동_검증결과_쉽게설명.txt'),
 'source_hashes_verified':len(manifest['files']),
 'measured_field_candidates':candidates,'matched_field_decisions':[{k:v for k,v in r.items() if k!='rows'} for r in matched],
 'weak_term_max_abs_delta_S_kB_per_spin':max(abs(c['delta_S_cold_kB_per_spin']) for c in weak['comparisons'])},ensure_ascii=False,indent=2))
