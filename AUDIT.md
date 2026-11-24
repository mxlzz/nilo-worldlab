# Auditoría de continuidad — 2026-11-24

Autoverificación documental del diseñador, no nueva reproducción científica. Origen fijado 0448ea6e5e02d77142cec270ea99ac0773493fc1; el corte dc1fdff166f75311ffcb7b72a3ac67c0f99cdb06 quedó superado por 12 commits. Revisión inicial: 42 alternativas, 18 ADR, cierres y resultados críticos; no revisión línea por línea de todos los archivos.

## Estado consolidado
E0: selección modular experimental (42 alternativas, seis familias ensayadas). E1: materia/owner, balances y calor/fase aislados. E2: dos leyes locales. E3: agregado granular causal H26/H28 contra nulo condicionado. E4: generación 2.5D. E5: suelo→contacto unilateral; H32 no pasó contacto aguas abajo. E6: reinicio de subsistemas y P0 bajo perfil limitado. Todos cierres provisionales. E7 abierta; E8 no iniciada.
Negativos: H19 nulo explica patrón; H21–23 estructura no validada; H36 soporte no robusto; H41 control de cruce inválido; H42 original no repetible, H42b geometría positiva contra nulo, no calor. P0-B nativa falló hashes por entorno, B-R aprobada con perfil común.
P0-C CI 45445405808; P0-D CI 45447643891: 4 casos, 125 granos, siete épocas, 12 estados de sufijo exactos Intel→AMD, consulta no mutante en caso probado, rechazo de petición duplicada. El workflow posterior 61866343 no tiene CI propia; evidencia corresponde a 9e1a6540. P0-E completo, GO condicional capacidad nueva y NO-GO motor final. Usuario aprobó P1.
P1-A original CI 45449986939 FAILURE por instrumento antes de FiPy. P1-A-R CI 45488606315 SUCCESS y WORLD_EVAL_JSON con gates vacíos, ver evidence/P1-A-R_READBACK.md. No afirmar energía integrada ni retroacción viva.

## Brechas
Tiempo/unidades nominales, falta de retroacción ambiental y diversidad funcional, portabilidad restringida, acople unilateral, 2.5D y ROI, rendimiento solo en escenarios pequeños, arneses específicos y artefactos temporales. Priorizar P1-B contrato material/energético. No crear scheduler, GUI, motor completo ni nuevos montículos sin evidencia. Chrono, BMI/HELICS, FiPy, modelos de agua o DEM/CFD solo si una necesidad demostrada lo justifica. Revisión heredada no exhaustiva.

La migración no valida automáticamente la física. Verificar cobertura, hashes, evidencia y exclusividad antes de habilitar ciencia. Si binarios no se recuperan, no afirmar preservación íntegra.

## Cierre WL-M
Ver [verificación](evidence/MIGRATION_VERIFICATION.md): 314/314 blobs heredados, 58/58 artefactos con digest, 487 ejecuciones inventariadas y 481 ZIP de logs. Seis fallos históricos sin jobs devuelven 404; se preservan metadatos y no bloquean P1-B.1. Simulación documental satisfactoria; no auditoría científica independiente. Siguiente tarea P1-B.1.
## P1-B.1 — autoverificación documental
Alineación: P1/E7, consecuencia ambiental atribuida a materia existente. Rechazada proyección gruesa como única solución espacial y lámina duplicada; propuesta partición fija de ROI con energía de columna completa. Flujos brutos nativos existen en fuente etiquetada Landlab 2.11.0; compatibilidad y balances no comprobados en ejecución. Conservación de entalpía no puede deducirse del delta masa. B.2a debe verificar campos, rutas, masas y export contra P0-D, sin modificar ley de erosión. Supuestos de mezcla/aislamiento y tiempos nominales declarados. GO documental únicamente, P1-B no aprobada; E7 abierta. Sin auditor independiente. Ver protocolo.

## P1-B.2a — GO instrumental
# P1-B.2a — Resultado y auditoría

CI [45502145469](https://github.com/mxlzz/nilo-worldlab/actions/runs/45502145469), job 138827906462, commit 41d52b3a00e16276c8f504d5d283292f8f22f2cb. SUCCESS y lectura WORLD_EVAL_JSON: cero gates fallidos, 72 pasos, 28/28 hashes de suelo exactos, K0 sin flujos, control H33 sin instrumentación exacto. Cuatro archivos NPZ con digests en JSON. Tiempo Python 1.384 s, job 40 s. Máximo residual nodal de masa 7.28e-10 kg; de flujo 5.03e-10 kg; rutas 4.66e-10 kg; residual export cero. Export total en casos activos 32.52 y 34.06 millones kg nominales; no equivale a energía transferida. Versiones y perfil correctos. Aceptación instrumental, no prueba de transporte energético.

Autoverificación: cumple P1-B.2a, sin modificación de física heredada. El protocolo no exige auditor externo para unidad. B.2b puede iniciarse solamente tras retención durable de NPZ y verificación de su digest. No usar tracker como verdad si no coincide con datos efectivos. Artefacto 14969549342, digest sha256:60028665102c3eb233b3c9d7ad15048be157e39d20189b69f831448bc76e635c, 782522 bytes, vence 2026-12-24T22:34:52Z. Se conserva JSON íntegro y se archiva ZIP/log en evidencia. No confundir success CI con éxito de P1-B.
