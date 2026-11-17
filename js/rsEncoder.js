/**
 * Codificación sistemática Reed-Solomon.
 *
 * Ver anteproyecto, sección 4.3: G(x) se construye con raíces alpha^1..alpha^(2t);
 * el mensaje se desplaza x^(2t) y se divide por G(x); el resto es la paridad.
 * T(x) = mensaje desplazado + resto -> los primeros k símbolos del codeword
 * son el mensaje original sin modificar (codificación sistemática).
 */

if (typeof require !== 'undefined') {
  var { polyMul, polyEval } = require('./poly.js');
}

/**
 * Construye el polinomio generador:
 *   G(x) = (x + alpha^1)(x + alpha^2)...(x + alpha^nsym)
 * (leading-first). nsym = 2t = cantidad de símbolos de paridad.
 */
function buildGeneratorPoly(gf, nsym) {
  let g = [1];
  for (let i = 1; i <= nsym; i++) {
    g = polyMul(gf, g, [1, gf.alphaTo(i)]);
  }
  return g;
}

/**
 * Codifica un mensaje de k símbolos (leading-first) con Reed-Solomon
 * sistemático, agregando nsym símbolos de paridad.
 *
 * Algoritmo: equivalente a "desplazar x^nsym · M(x) y dividir por G(x)",
 * implementado como una división polinómica en un solo recorrido (idéntico
 * en espíritu al registro de desplazamiento usado para CRC).
 *
 * @param {GaloisField} gf
 * @param {number[]} msgIn - símbolos del mensaje, leading-first, largo k
 * @param {number} nsym - cantidad de símbolos de paridad (2t)
 * @param {number[]} [generator] - generador ya calculado (opcional, para reutilizar)
 * @returns {{codeword:number[], parity:number[], generator:number[]}}
 */
function rsEncode(gf, msgIn, nsym, generator) {
  if (!Array.isArray(msgIn) || msgIn.length === 0) {
    throw new Error('El mensaje debe ser un array no vacío de símbolos');
  }
  const gen = generator || buildGeneratorPoly(gf, nsym);

  // msgOut = mensaje desplazado x^nsym (se agregan nsym ceros al final)
  const msgOut = msgIn.concat(new Array(nsym).fill(0));

  // División polinómica en un solo recorrido: en cada posición i, se
  // "cancela" el coeficiente líder actual restando (XOR) coef * G(x)
  // desplazado a esa posición. Como G(x) es mónico (gen[0] = 1), la
  // posición i queda en 0, exactamente el efecto de una división larga.
  for (let i = 0; i < msgIn.length; i++) {
    const coef = msgOut[i];
    if (coef !== 0) {
      for (let j = 0; j < gen.length; j++) {
        msgOut[i + j] ^= gf.mul(gen[j], coef);
      }
    }
  }

  // Lo que queda en las últimas nsym posiciones es el resto = paridad.
  const parity = msgOut.slice(msgIn.length);

  return { codeword: msgIn.concat(parity), parity, generator: gen };
}

/**
 * Verifica un codeword evaluándolo en las raíces del generador
 * (alpha^1..alpha^nsym): si el código es válido, todas las evaluaciones
 * dan 0 (síndromes nulos). Útil como chequeo rápido tras codificar.
 */
function rsVerifyCodeword(gf, codeword, nsym) {
  for (let i = 1; i <= nsym; i++) {
    if (polyEval(gf, codeword, gf.alphaTo(i)) !== 0) return false;
  }
  return true;
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { buildGeneratorPoly, rsEncode, rsVerifyCodeword };
}
