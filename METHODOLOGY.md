# Metodología experimental

Aplicar inherited/WORLD_RESEARCH/PROCEDIMIENTO_Y_AUDITORIA.md y los criterios E0–E8; no repetir 42 alternativas salvo brecha nueva. Antes de implementar, contrastar soluciones existentes, opción mínima/no implementar, evidencia original, licencias, dependencias y costo.

## Protocolo por unidad
- Pregunta vinculada al objetivo/hito, hipótesis falsable y nula, supuestos.
- Prerregistrar semillas, versiones, parámetros, entradas, controles (negativo, ablación, sham, plantilla o nulo según pregunta), métrica, umbral, presupuesto y parada. Si controles no discriminan, resultado inconcluso/no aprobado.
- Ejecutar mínima prueba y conservar código, entorno, comandos, logs, datos, resultados negativos y fallos instrumentales. No ajustar gates tras conocer resultados. CI SUCCESS no basta: cotejar criterios.
- Verificar masa, energía, identidad, geometría, causalidad, tiempo y unidades según corresponda. Distinguir unidad, laboratorio aislado, composición integrada, generalización y física calibrada. No inferir p iid de MCMC condicionado.
- Auditar valor marginal y alineación. Elegir continuar/corregir/replantear/retroceder/detener; preferir auditor distinto, declarar autoverificación si coincide.

## Contratos heredados
WorldEngine pin 3b9113ecfad884cab4b32081d138a15e6262ada8, alturas arbitrarias convertidas a nominales. Landlab 2.11 y PyBullet 3.2.7 separados, xarray propietario lógico, NetCDF y cache Bullet. No duplicar stocks. Terreno→contacto unilateral; no asumir retroacción. Relojes geomorfológico benchmark, Bullet y térmico independientes; no renombrar SI. Perfil P0-B-R: Python 3.13.18/3.11.18, NumPy 2.5.5/2.4.6, xarray 2026.9.0, SciPy 1.18.1/1.17.1; NPY_DISABLE_CPU_FEATURES=X86_V3,X86_V4,AVX512_ICL,AVX512_SPR y OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1. Exactitud solo bajo perfiles/hosts ensayados; no redondear referencias históricas.

## Persistencia
Probar restart completo sin replay del pasado, corrupción, identidad, cache y ausencia de futuro entregado al consumidor. Consultas no mutantes, replay sin instrumentación. Artefactos temporales: archivo con digest o regeneración comprobada, no URL eterna. Si faltan evidencias críticas, bloquear unidad afectada. No gasto pago sin autorización. Leer HEAD antes de cambios, resolver conflictos, verificar escrituras, actualizar STATE.json al final. No almacenar secretos.