# mPASS 알고리즘 현대 C++ 재작성 분석

## 분석 대상

`original_source/cal.cpp` 및 관련 파일 (`pass.h`, `cal.h`, `reader.cpp`, `soga.cpp`, `pass.cpp`)

---

## Step 1 – 로직 분석 및 의사 코드(Pseudo-code)

### 전체 파이프라인 요약

mPASS(modified PASS)는 단백질 결합 부위(Binding Site)를 탐지하기 위한 알고리즘입니다.
단백질 표면에 구 형태의 "프로브(Probe)"를 빈 공간에 점점 채워넣고,
고밀도로 모인 프로브 군집(Cluster)을 결합 부위 후보로 분류합니다.

```
FUNCTION mPASS(pdb_file):

  # 1. 파일 로딩
  atomProps  ← readAtomPropertyFile("atom_property")
  gridProps  ← readGridPropertyFile("grid_property")
  protein    ← readPdbFile(pdb_file, atomProps)      // ATOM/HETATM 레코드 파싱
  
  # 2. 공간 분할 (Cell list)
  grid ← buildCellList(protein, cellSize=2.0Å)       // 3D 균일 격자 구축
  
  # 3. 1번째 레이어 프로브 생성
  FOR EACH triple (i, j, k) of near protein atoms:
    FOR EACH candidate position p IN getTripleVertex(i, j, k, probeRadius):
      IF bumpCheckAndScore(p, protein, grid) PASSES:
        ADD new Probe(p) to probeList
  
  # 4. 초기 DBSCAN 클러스터링 + 겹침 제거
  clusterProbes(probeList, eps=0.7, minPts=1)
  weedOverlappingProbes(probeList, weedRadius=0.7)    // 너무 가까운 프로브 제거
  
  # 5. BC(Buriedness Count) 계산 & 임계값 도출
  FOR EACH surviving probe p:
    p.buriedness ← countAtomsWithin14A(p, grid)
  deriveBcCutoff(probeList, proteinSize)
  
  # 6. 2~4번째 레이어 프로브 생성 (각 레이어는 이전 레이어 + 단백질 원자 사용)
  FOR layer = 2 TO 4:
    FOR EACH triple (a, b, c) using probes from previous layers:
      FOR EACH candidate position p IN getTripleVertex(a, b, c, probeRadius):
        IF bumpCheckAndScore(p, protein, grid) PASSES:
          ADD new Probe(p) to probeList
  
  # 7. 최종 DBSCAN 클러스터링
  assignProbeNeighbors(probeList, cutoff=2.8Å)
  clusterProbes(probeList, eps=2.5, minPts=16)
  buildClusterSet(probeList) → clusterMap
  
  # 8. 클러스터 특성 계산
  FOR EACH cluster c:
    c.depth   ← minDist(opening probes → deepest probe)
    c.width   ← maxPairwiseDist(all probes in c)
    c.opening ← dbscanClusteringForSurface(topProbes, eps=4, minPts=1)
  
  # 9. PLB(결합 가능성) 점수 계산
  FOR EACH cluster c:
    c.plb   ← SUM(sogaIndex  for residue in c.contacts)
    c.hydro ← SUM(kyte-doolittle for residue in c.contacts)
  
  RETURN clusterMap
```

---

### 핵심 알고리즘 1 – GetTripleVertex (3구 교점 계산)

세 단백질 원자를 동시에 접하는 프로브 위치를 구합니다.

```
FUNCTION getTripleVertex(centerI, centerJ, centerK, rI, rJ, rK, probeR):

  # 로컬 좌표계 구성
  Xu ← normalize(centerJ − centerI)        // X축
  Zu ← normalize(Xu × (centerK − centerI)) // 평면의 법벡터 (Z축)
  Yu ← Zu × Xu                             // Y축

  # 이분 투영점 계산 (각 변의 수직이등분 평면)
  t_ij ← (‖rI+rp‖² − ‖rJ+rp‖²) / (2 × ‖centerJ−centerI‖²)
  T_ij ← midpoint(centerI, centerJ) + (centerJ−centerI) × t_ij

  t_ik ← (‖rI+rp‖² − ‖rK+rp‖²) / (2 × ‖centerK−centerI‖²)
  T_ik ← midpoint(centerI, centerK) + (centerK−centerI) × t_ik

  # 평면 내 기저점 Rb 계산
  t   ← dot(T_ik−T_ij, T_ik−centerI) / dot(T_ik−centerI, Yu)
  Rb  ← T_ij + Yu × t

  # 평면 외 높이 h 계산
  hSq ← (rI+rp)² − ‖Rb − centerI‖²
  IF hSq < 0: RETURN null  // 실수해 없음

  RETURN [Rb + Zu×√hSq,   // 위쪽 해
          Rb − Zu×√hSq]   // 아래쪽 해
```

---

### 핵심 알고리즘 2 – DBSCAN 클러스터링

```
FUNCTION clusterProbes(probeList, eps, minPts, startId):

  # Step 1: 밀도 계산
  FOR EACH probe i IN probeList:
    i.density ← 0; i.mergedProbes ← []
    FOR EACH neighbor j IN i.nearProbes (미리 계산됨):
      IF j IN probeList AND dist(i,j) < eps:
        i.density++
        i.mergedProbes.append(j)

  # Step 2: 클러스터 확장
  clusterId ← startId
  FOR EACH unvisited probe i IN probeList:
    MARK i as visited
    IF i.density < minPts: CONTINUE  // 노이즈

    clusterId++; i.clusterId ← clusterId
    queue ← i.mergedProbes

    FOR EACH k IN queue:
      IF k is unvisited:
        MARK k as visited
        IF k.density >= minPts:
          ADD k.mergedProbes TO queue  // BFS 확장
      IF k.clusterId == 0:
        k.clusterId ← clusterId

  RETURN clusterId
```

---

### 핵심 알고리즘 3 – BC(Buriedness Count) 계산

```
FUNCTION calculateBuriedness(probePos, grid):
  wholesaleCount ← 0
  exactCount     ← 0

  FOR EACH neighboring cell c IN grid:
    IF c.maxDist ≤ 14.0:           // 셀 전체가 반경 내에 있음
      wholesaleCount += c.atoms.count
    ELSE:
      FOR EACH atom a IN c:
        IF dist(probePos, a) ≤ 14.0:
          exactCount++

  RETURN wholesaleCount + exactCount
```

---

## Step 2 – 테스트 설계

| # | 테스트명 | 범주 | 입력 | 기대 출력 |
|---|---|---|---|---|
| 1 | Vector3 산술 | Happy | a=(1,2,3), b=(4,5,6) | sum=(5,7,9), dot=32, cross=(-3,6,-3) |
| 2 | 유클리드 거리 | Happy | (0,0,0), (3,4,0) | 5.0 |
| 3 | 정규화 | Happy | (0,0,5) | (0,0,1), length=1 |
| 4 | 영벡터 정규화 | Edge | (0,0,0) | domain_error 예외 |
| 5 | 직교 벡터 각도 | Happy | (+x, +y) | 90° |
| 6 | PDB 파싱 – 기본 | Happy | ALA 1개 + GLY 1개 | 3원자, 2잔기 |
| 7 | PDB 파싱 – 바운딩박스 | Happy | x=[1,5] | minX=1-5, maxX=5+5 |
| 8 | PDB 파싱 – 미지의 원자 | Edge | CB 없는 맵 | CB 줄 스킵됨 |
| 9 | PDB 파싱 – 빈 파일 | Edge | 빈 스트림 | 원자 0개, 잔기 0개 |
| 10 | GetTripleVertex – 정삼각형 | Happy | 3개 동일구(r=1.7, rp=0.7) | 거리 오차 ≤ 0.02 Å |
| 11 | GetTripleVertex – 동일선상 | Edge | 3점 일직선 | nullopt |
| 12 | GetTripleVertex – 동일점 | Edge | A=B | nullopt |
| 13 | DBSCAN – 5개 밀집 프로브 | Happy | 간격 0.5Å, eps=1.5 | 클러스터 1개 |
| 14 | DBSCAN – 2개 분리된 그룹 | Happy | 간격 0.5Å, 먼 거리 | 클러스터 2개 |
| 15 | DBSCAN – 단독 프로브 (minPts=2) | Edge | 프로브 1개 | 노이즈(cid=0) |
| 16 | DBSCAN – 빈 입력 | Edge | 빈 벡터, startId=42 | 42 반환 |
| 17 | BC 임계값 – 자동 계산 | Happy | 500원자, top10=200 | bcCutoff ∈ [100, 160] |
| 18 | BC 임계값 – Override | Happy | override=0.65 | bcCutoff=65 |
| 19 | SOGA 테이블 – 20종 로딩 | Happy | init 호출 | TRP=2.518, ILE=4.5 |
| 20 | PLB 점수 – ALA+ILE | Happy | residuesContact=[ALA,ILE] | plb=1.707, hydro=6.3 |
| 21 | PLB 점수 – 미지 잔기 | Edge | residuesContact=[UNK] | plb=0, hydro=0 |

---

## Step 3 – 현대 C++17 재작성 결과물

| 파일 | 역할 | 주요 개선 사항 |
|---|---|---|
| `geometry.hpp` | Vector3 수학 | `[[nodiscard]]`, `constexpr`, `std::clamp`, 예외 처리 |
| `data_types.hpp` | 데이터 구조 | `enum class Polarity`, `constexpr` 상수, 스마트 포인터 없이 값 타입 |
| `grid.hpp` | 공간 격자 | 함수로 캡슐화, 단일 코드 경로 |
| `bc_calculator.hpp` | BC 계산 | `std::clamp`, 람다 없는 명확한 로직 |
| `probe_placer.hpp` | 프로브 배치 | `std::optional`, 영벡터 검사 |
| `dbscan.hpp` | DBSCAN 클러스터링 | 템플릿 BFS, `buildClusterSet` 분리 |
| `residue_properties.hpp` | SOGA/소수성 | 정적 초기화 테이블 |
| `file_reader.hpp` | 파일 파서 | `std::stod` 대신 `atof`, `std::istream&` 인터페이스 |
| `pass_modern.cpp` | 메인 진입점 | 파이프라인 명시적 오케스트레이션 |
| `tests_modern.cpp` | 단위 테스트 | 23개 테스트, 프레임워크 없이 자체 실행 |

---

## Step 4 – 검증 가이드 및 언어 특성 차이

### 4-1. 논리 동일성 확인 방법

```bash
# 빌드
cd modern_source
g++ -std=c++17 -O2 -I. tests_modern.cpp -o run_tests && ./run_tests

# 기존 코드와 같은 PDB 파일로 출력 비교 (probe 좌표, cluster 수, BC 값)
diff <(./original_pass input.pdb 2>&1) <(./pass_modern input.pdb 2>&1)
```

### 4-2. 언어 특성상 차이점

| 항목 | 원본 코드 | 현대 코드 | 영향 |
|---|---|---|---|
| **정수 오버플로우** | `int` 로 BC 누적 (최대 ~수만) | 동일 `int` — 오버플로우 없음 | 없음 |
| **부동소수점 반올림** | `double` 산술 | `double` 산술 동일 | 없음 |
| **atof vs std::stod** | `atof("")` = 0.0 (UB 없음) | `std::stod("")` → 예외 포착 후 0.0 | 행동 동일, 오류 감지 향상 |
| **전역 변수** | `g_probes`, `g_clusterMap` 등 전역 | `SimulationState` 구조체로 묶음 | 스레드 안전성 향상 |
| **원시 포인터** | `new GridPoint`, `new Probe` | 값 타입 `std::vector<Probe>` | 메모리 누수 제거 |
| **매크로 상수** | `#define PI 3.14159...` | `constexpr double PI = ...` | 타입 안전성 |
| **SQUARE 매크로** | `#define SQUARE(x) ((x)*(x))` | `template sq(x)` | 타입 안전성, 인라인 |
| **정렬 순서** | `std::sort` with raw comparator | `std::sort` with lambda | 동일 |
| **map vs unordered_map** | `std::map<int, Cluster*>` | `std::map<int, Cluster>` | 값 의미론, 소유권 명확 |

### 4-3. 테스트로 확인해야 할 수치적 동일성

1. **GetTripleVertex**: `|dist(probe, atom) - (r_atom + r_probe)| ≤ 0.01 Å` (원본과 동일 허용 오차)
2. **BC 값**: 동일한 PDB + 동일한 grid_property 파일 → 동일한 정수 BC 값
3. **bcCutoff**: `top10Mean * ratio` – 원본과 1e-6 이내
4. **PLB 점수**: 잔기 목록이 동일하면 소수점 4자리까지 일치
5. **클러스터 수**: 동일한 DBSCAN 파라미터(eps=2.5, minPts=16) → 동일한 클러스터 수

### 4-4. 알려진 버그 수정 (원본 대비)

`weedFirstLayer()` 내 `density` 비교 버그:

```cpp
// 원본 (버그): a와 a를 비교
if(g_probes[a]->density > g_probes[a]->density) { ... }

// 수정된 코드: a와 b를 비교
if (state.probes[i].density > other.density) { ... }
```

이 버그로 인해 원본에서는 `density` 비교 분기가 항상 `else` 경로로 진입합니다.
현대 코드는 올바른 비교를 수행합니다.
