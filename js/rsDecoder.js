/**
 * Decodificación Reed-Solomon: síndromes, polinomio localizador de errores
 * (Berlekamp-Massey), localización (búsqueda de Chien), magnitud del error
 * (algoritmo de Forney) y corrección.
 *
 * Ver anteproyecto, sección 4.4. Convenciones usadas en este archivo:
 *  - `received` (codeword recibido) es leading-first, igual que en rsEncoder.js:
 *    índice 0 = coeficiente de x^(n-1), índice n-1 = coeficiente de x^0.
 *  - Los síndromes se calculan en las raíces alpha^1 .. alpha^(2t) del
 *    generador (convención "b=1"), como en el anteproyecto.
 *  - El polinomio localizador Λ(x) y el evaluador de errores Ω(x) se manejan
 *    en orden ASCENDENTE (índice 0 = coeficiente de x^0), que es la
 *    convención habitual del algoritmo de Berlekamp-Massey.
 *  - "posición" de un error se identifica por su EXPONENTE de x (0..n-1),
 *    igual que en el anteproyecto ("error en la posición x^3"). Cada
 *    posición reportada incluye también el índice de array equivalente
 *    (index = n-1-exponente) para aplicar la corrección directamente sobre
 *    el array `received`.
 */

if (typeof require !== 'undefined') {
  var { polyMul, polyEval, polyEvalAscending } = require('./poly.js');
}

/**
 * Síndromes S_1..S_nsym: S[i-1] = received(alpha^i), evaluado con Horner
 * sobre el codeword leading-first.
 */
function computeSyndromes(gf, received, nsym) {
  const syndromes = new Array(nsym);
  for (let i = 1; i <= nsym; i++) {
    syndromes[i - 1] = polyEval(gf, received, gf.alphaTo(i));
  }
  return syndromes;
}

/** true si todos los síndromes son 0 (codeword sin errores detectables). */
function syndromesAllZero(syndromes) {
  return syndromes.every((s) => s === 0);
}

/**
 * Algoritmo de Berlekamp-Massey: a partir de los síndromes, encuentra el
 * polinomio localizador de errores Λ(x) de menor grado que satisface la
 * ecuación clave. Devuelve Λ(x) en orden ascendente (Λ[0] = 1 siempre) y la
 * cantidad de errores detectados (= grado de Λ(x) = longitud del LFSR).
 *
 * Es la versión general y eficiente del método usado a mano en el
 * anteproyecto con Peterson-Gorenstein-Zierler (PGZ): para pocos errores
 * ambos dan el mismo resultado, pero Berlekamp-Massey no requiere adivinar
 * de antemano cuántos errores hay ni resolver sistemas de ecuaciones cada
 * vez, lo que lo hace preferible para la implementación en la calculadora.
 */
function berlekampMassey(gf, syndromes) {
  const n = syndromes.length; // = nsym = 2t

  let C = [1]; // Λ(x) candidato actual (ascendente)
  let B = [1]; // Λ(x) candidato de la iteración previa al último cambio de longitud
  let L = 0; // cantidad de errores detectados hasta el momento
  let m = 1; // pasos desde el último cambio de longitud
  let b = 1; // discrepancia en el último cambio de longitud

  const ensureLength = (arr, len) => {
    while (arr.length < len) arr.push(0);
  };

  for (let nIter = 0; nIter < n; nIter++) {
    // discrepancia: delta = S[nIter] + sum_{i=1}^{L} C[i] * S[nIter-i]
    let delta = syndromes[nIter];
    for (let i = 1; i <= L; i++) {
      delta ^= gf.mul(C[i] || 0, syndromes[nIter - i]);
    }

    if (delta === 0) {
      m += 1;
      continue;
    }

    const coef = gf.div(delta, b);

    if (2 * L <= nIter) {
      const T = C.slice(); // copia de C antes de modificarlo
      ensureLength(C, m + B.length);
      for (let i = 0; i < B.length; i++) {
        C[i + m] ^= gf.mul(coef, B[i]);
      }
      L = nIter + 1 - L;
      B = T;
      b = delta;
      m = 1;
    } else {
      ensureLength(C, m + B.length);
      for (let i = 0; i < B.length; i++) {
        C[i + m] ^= gf.mul(coef, B[i]);
      }
      m += 1;
    }
  }

  return { locator: C, numErrors: L };
}

/**
 * Búsqueda de Chien: prueba, para cada posible posición 0..n-1 del
 * codeword, si alpha^{-exponente} es raíz de Λ(x). Donde lo es, hay un
 * error en esa posición.
 */
function chienSearch(gf, locator, n) {
  const positions = [];
  for (let exponent = 0; exponent < n; exponent++) {
    const X = gf.alphaTo(exponent); // valor localizador X_l = alpha^exponente
    const root = gf.inv(X); // raíz que se busca: X_l^{-1}
    if (polyEvalAscending(gf, locator, root) === 0) {
      positions.push({
        exponent, // "posición x^exponente" (como en el anteproyecto)
        index: n - 1 - exponent, // índice equivalente en el array leading-first
        X,
        root,
      });
    }
  }
  return positions;
}

/**
 * Derivada formal de un polinomio ascendente sobre un cuerpo de
 * característica 2: los términos de exponente par se anulan, y los de
 * exponente impar bajan un grado con el mismo coeficiente (no se multiplica
 * por el exponente: en GF(2^m), multiplicar por un entero impar equivale a
 * dejar el coeficiente igual).
 */
function formalDerivative(coeffsAscending) {
  const result = new Array(Math.max(coeffsAscending.length - 1, 0)).fill(0);
  for (let i = 1; i < coeffsAscending.length; i += 2) {
    result[i - 1] = coeffsAscending[i];
  }
  return result;
}

/**
 * Algoritmo de Forney: calcula la magnitud de cada error ya localizado.
 * Ω(x) = [S(x)·Λ(x)] mod x^nsym (evaluador de errores); con raíces
 * alpha^1..alpha^nsym (convención b=1), la magnitud es
 * Y_l = Ω(X_l^{-1}) / Λ'(X_l^{-1}) (sin signo negativo: en GF(2^m), -a = a).
 */
function forneyErrorValues(gf, syndromes, locator, errorPositions, nsym) {
  const lambdaPrime = formalDerivative(locator);
  const fullProduct = polyMul(gf, syndromes, locator); // syndromes y locator, ambos ascendentes
  const omega = fullProduct.slice(0, nsym); // mod x^nsym: nos quedamos con los términos de grado < nsym

  return errorPositions.map(({ root, exponent, index }) => {
    const omegaAtRoot = polyEvalAscending(gf, omega, root);
    const lambdaPrimeAtRoot = polyEvalAscending(gf, lambdaPrime, root);
    if (lambdaPrimeAtRoot === 0) {
      throw new Error(
        "Λ'(x) se anula en una raíz del localizador: caso degenerado, no se puede aplicar Forney"
      );
    }
    const magnitude = gf.div(omegaAtRoot, lambdaPrimeAtRoot);
    return { index, exponent, magnitude };
  });
}

/**
 * Decodifica un codeword recibido: detecta, localiza y corrige errores.
 *
 * @param {GaloisField} gf
 * @param {number[]} received - codeword recibido, leading-first, largo n
 * @param {number} nsym - símbolos de paridad del código (2t)
 * @returns {{
 *   success: boolean,
 *   hadErrors: boolean,
 *   corrected: number[],
 *   syndromes: number[],
 *   locator: number[],
 *   numErrorsFound: number,
 *   errors: {index:number, exponent:number, magnitude:number}[],
 *   reason?: string,
 * }}
 */
function rsDecode(gf, received, nsym) {
  const n = received.length;
  const t = Math.floor(nsym / 2);
  const syndromes = computeSyndromes(gf, received, nsym);

  if (syndromesAllZero(syndromes)) {
    return {
      success: true,
      hadErrors: false,
      corrected: received.slice(),
      syndromes,
      locator: [1],
      numErrorsFound: 0,
      errors: [],
    };
  }

  const { locator, numErrors } = berlekampMassey(gf, syndromes);

  if (numErrors === 0 || numErrors > t) {
    return {
      success: false,
      hadErrors: true,
      corrected: received.slice(),
      syndromes,
      locator,
      numErrorsFound: numErrors,
      errors: [],
      reason:
        numErrors > t
          ? `Se detectaron más errores (${numErrors}) que la capacidad de corrección del código (t=${t}). Falla de decodificación: el codeword no es corregible.`
          : 'Los síndromes son inconsistentes con un patrón de error válido (Λ(x) trivial pero síndromes no nulos).',
    };
  }

  const errorPositions = chienSearch(gf, locator, n);

  if (errorPositions.length !== numErrors) {
    return {
      success: false,
      hadErrors: true,
      corrected: received.slice(),
      syndromes,
      locator,
      numErrorsFound: numErrors,
      errors: [],
      reason:
        'La búsqueda de Chien no encontró raíces válidas para todas las posiciones de error: ' +
        'patrón no corregible (probablemente excede t o el receptor recibió un codeword inválido).',
    };
  }

  const errors = forneyErrorValues(gf, syndromes, locator, errorPositions, nsym);

  const corrected = received.slice();
  for (const err of errors) {
    corrected[err.index] ^= err.magnitude;
  }

  // Verificación final: el codeword corregido debe tener síndromes nulos.
  const verifySyndromes = computeSyndromes(gf, corrected, nsym);
  const verified = syndromesAllZero(verifySyndromes);

  return {
    success: verified,
    hadErrors: true,
    corrected,
    syndromes,
    locator,
    numErrorsFound: numErrors,
    errors,
    reason: verified
      ? undefined
      : 'La corrección no pasó la verificación final de síndromes (patrón de error no corregible con este código).',
  };
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    computeSyndromes,
    syndromesAllZero,
    berlekampMassey,
    chienSearch,
    formalDerivative,
    forneyErrorValues,
    rsDecode,
  };
}
