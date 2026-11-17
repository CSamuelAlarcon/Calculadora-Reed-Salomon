/**
 * Operaciones con polinomios sobre un cuerpo GF(2^m).
 *
 * Convención: los polinomios se representan como arrays de coeficientes en
 * orden "leading-first" (índice 0 = coeficiente del término de mayor grado),
 * igual que se escribieron a mano en el anteproyecto — por ejemplo,
 * T(x) = [12, 6, 3, 3, 12, 8, 0] para los coeficientes de x^6 a x^0.
 */

/** Grado del polinomio (índice del último coeficiente). */
function polyDegree(p) {
  return p.length - 1;
}

/** Elimina ceros líderes sobrantes (normaliza el polinomio a su grado real). */
function polyTrim(p) {
  let i = 0;
  while (i < p.length - 1 && p[i] === 0) i++;
  return p.slice(i);
}

/** Suma/resta de polinomios (en GF(2^m) son la misma operación: XOR). */
function polyAdd(p, q) {
  const len = Math.max(p.length, q.length);
  const result = new Array(len).fill(0);
  for (let i = 0; i < p.length; i++) result[i + (len - p.length)] ^= p[i];
  for (let i = 0; i < q.length; i++) result[i + (len - q.length)] ^= q[i];
  return result;
}

/** Multiplicación de polinomios: convolución de coeficientes sobre GF(2^m). */
function polyMul(gf, p, q) {
  const result = new Array(p.length + q.length - 1).fill(0);
  for (let j = 0; j < q.length; j++) {
    if (q[j] === 0) continue;
    for (let i = 0; i < p.length; i++) {
      if (p[i] === 0) continue;
      result[i + j] ^= gf.mul(p[i], q[j]);
    }
  }
  return result;
}

/** Escala todos los coeficientes de p por una constante k. */
function polyScale(gf, p, k) {
  return p.map((c) => gf.mul(c, k));
}

/** Evalúa el polinomio p (leading-first) en x, por el método de Horner. */
function polyEval(gf, p, x) {
  let result = p[0];
  for (let i = 1; i < p.length; i++) {
    result = gf.mul(result, x) ^ p[i];
  }
  return result;
}

/**
 * Evalúa un polinomio en convención ASCENDENTE (índice 0 = coeficiente de
 * x^0, es decir p[i] es el coeficiente de x^i). Se usa para el localizador
 * de errores Λ(x) y el evaluador de errores Ω(x) del decodificador, que por
 * convención de Berlekamp-Massey se manejan en este orden. Reutiliza
 * polyEval invirtiendo el array (mismo resultado, sin duplicar lógica).
 */
function polyEvalAscending(gf, p, x) {
  return polyEval(gf, p.slice().reverse(), x);
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    polyDegree,
    polyTrim,
    polyAdd,
    polyMul,
    polyScale,
    polyEval,
    polyEvalAscending,
  };
}
