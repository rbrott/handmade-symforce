"use strict";

/* global document, window, twgl, m3 */

const canvas = document.querySelector("canvas");
const gl = canvas.getContext("webgl");

const vs = `
attribute vec2 a_position;
uniform mat3 u_matrix;
void main() {
  gl_Position = vec4((u_matrix * vec3(a_position, 1)).xy, 0, 1);
}
`;

const fs = `
precision mediump float;
uniform vec4 u_color;
void main() {
  gl_FragColor = u_color;
}
`;

// compiles shaders, links program, looks up locations
const programInfo = twgl.createProgramInfo(gl, [vs, fs]);

// calls gl.createBuffer, gl.bindBuffer, gl.bufferData
let bufferInfo;

const camera = {
  x: 0,
  y: 0,
  zoom: 1
};

let viewProjectionMat;

function makeCameraMatrix() {
  const zoomScale = 1 / camera.zoom;
  let cameraMat = m3.identity();
  cameraMat = m3.translate(cameraMat, camera.x, camera.y);
  cameraMat = m3.scale(cameraMat, zoomScale, zoomScale);
  return cameraMat;
}

function updateViewProjection() {
  // same as ortho(0, width, height, 0, -1, 1)
  const projectionMat = m3.projection(gl.canvas.width, gl.canvas.height);
  const cameraMat = makeCameraMatrix();
  let viewMat = m3.inverse(cameraMat);
  viewProjectionMat = m3.multiply(projectionMat, viewMat);
}

function draw() {
  gl.clear(gl.COLOR_BUFFER_BIT);

  updateViewProjection();

  gl.useProgram(programInfo.program);

  const x = 0;
  const y = 0;
  const scale = 1;
  const color = [0, 0, 1, 1];

  // calls gl.bindBuffer, gl.enableVertexAttribArray, gl.vertexAttribPointer
  twgl.setBuffersAndAttributes(gl, programInfo, bufferInfo);

  let mat = m3.identity();
  mat = m3.translate(mat, x, y);
  mat = m3.scale(mat, scale, scale);

  // calls gl.uniformXXX
  twgl.setUniforms(programInfo, {
    u_matrix: m3.multiply(viewProjectionMat, mat),
    u_color: color
  });

  // calls gl.drawArrays or gl.drawElements
  twgl.drawBufferInfo(gl, bufferInfo);
}

function getClipSpaceMousePosition(e) {
  // get canvas relative css position
  const rect = canvas.getBoundingClientRect();
  const cssX = e.clientX - rect.left;
  const cssY = e.clientY - rect.top;

  // get normalized 0 to 1 position across and down canvas
  const normalizedX = cssX / canvas.clientWidth;
  const normalizedY = cssY / canvas.clientHeight;

  // convert to clip space
  const clipX = normalizedX * 2 - 1;
  const clipY = normalizedY * -2 + 1;

  return [clipX, clipY];
}

let startInvViewProjMat;
let startCamera;
let startPos;
let startClipPos;
let startMousePos;

function moveCamera(e) {
  const pos = m3.transformPoint(
    startInvViewProjMat,
    getClipSpaceMousePosition(e)
  );

  camera.x = startCamera.x + startPos[0] - pos[0];
  camera.y = startCamera.y + startPos[1] - pos[1];
  draw();
}

function handleMouseMove(e) {
  moveCamera(e);
}

function handleMouseUp(e) {
  draw();
  window.removeEventListener("mousemove", handleMouseMove);
  window.removeEventListener("mouseup", handleMouseUp);
}

canvas.addEventListener("mousedown", (e) => {
  e.preventDefault();
  window.addEventListener("mousemove", handleMouseMove);
  window.addEventListener("mouseup", handleMouseUp);

  startInvViewProjMat = m3.inverse(viewProjectionMat);
  startCamera = Object.assign({}, camera);
  startClipPos = getClipSpaceMousePosition(e);
  startPos = m3.transformPoint(startInvViewProjMat, startClipPos);
  startMousePos = [e.clientX, e.clientY];
  draw();
});

canvas.addEventListener("wheel", (e) => {
  e.preventDefault();
  const [clipX, clipY] = getClipSpaceMousePosition(e);

  // position before zooming
  const [preZoomX, preZoomY] = m3.transformPoint(
    m3.inverse(viewProjectionMat),
    [clipX, clipY]
  );

  // multiply the wheel movement by the current zoom level
  // so we zoom less when zoomed in and more when zoomed out
  const newZoom = camera.zoom * Math.pow(2, e.deltaY * -0.01);
  camera.zoom = Math.max(0.02, Math.min(100, newZoom));

  updateViewProjection();

  // position after zooming
  const [postZoomX, postZoomY] = m3.transformPoint(
    m3.inverse(viewProjectionMat),
    [clipX, clipY]
  );

  // camera needs to be moved the difference of before and after
  camera.x += preZoomX - postZoomX;
  camera.y += preZoomY - postZoomY;

  draw();
});

fetch("sparsity.bin").then((response) => {
  return response.arrayBuffer();
}).then((buffer) => {
  // i32 - dim
  // i32 - nnz
  // i32[dim + 1] - col_starts
  // i32[nnz] - row_indices
  const dim = new Int32Array(buffer, 0, 1)[0];
  const nnz = new Int32Array(buffer, 4, 1)[0];
  const colStarts = new Int32Array(buffer, 8, dim + 1);
  const rowIndices = new Int32Array(buffer, 8 + (dim + 1) * 4, nnz);

  const data = new Float32Array(8 * nnz);
  const indices = new Uint32Array(6 * nnz);

  gl.getExtension("OES_element_index_uint");

  let nz_index = 0;
  for (let col = 0; col < dim; ++col) {
    while (nz_index < colStarts[col + 1]) {
      const row = rowIndices[nz_index];
      data[8 * nz_index] = col;
      data[8 * nz_index + 1] = row;
      data[8 * nz_index + 2] = col + 1;
      data[8 * nz_index + 3] = row;
      data[8 * nz_index + 4] = col + 1;
      data[8 * nz_index + 5] = row + 1;
      data[8 * nz_index + 6] = col;
      data[8 * nz_index + 7] = row + 1;

      indices[6 * nz_index] = 4 * nz_index;
      indices[6 * nz_index + 1] = 4 * nz_index + 1;
      indices[6 * nz_index + 2] = 4 * nz_index + 2;
      indices[6 * nz_index + 3] = 4 * nz_index;
      indices[6 * nz_index + 4] = 4 * nz_index + 2;
      indices[6 * nz_index + 5] = 4 * nz_index + 3;

      nz_index++;
    }
  }

  bufferInfo = twgl.createBufferInfoFromArrays(gl, {
    a_position: {
      numComponents: 2,
      data,
    },
    indices,
  });

  draw();
});
