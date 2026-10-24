# The two films

* **The science and the steps** (`explainer`): what TripleSpec records, what the atmosphere and the
  sky do to it, and what each step of the reduction does about it, with the program's own data.
* **Hands-on tutorial** (`tutorial`): NGC 134 reduced from the raw frames in the interface, as a
  narrated screen recording; long computations are shown faster, with a badge saying how much.

Each is made in the same way:

| step | script | writes |
|---|---|---|
| narration | `voice.py <film>` (Chatterbox Turbo, a cloned voice; word times by forced alignment) | `build/<film>/voice/` |
| the recording (tutorial) | `record_gui.py` (drives the real interface off screen) | `build/tutorial/rec/` |
| the data (explainer) | `export_data.py` | `film/public/explainer/` |
| one clock | `assemble_explainer.py`, `assemble_tutorial.py` | `film/public/<film>/timeline.json`, `build/<film>/placements.json` |
| the soundtrack | `audio.py <film>` | `build/<film>/mix.wav` |
| the picture | `serve.sh` (Vite, port 5241), then `render.py <film>` | `build/<film>/TripleSpec_Reducer_<film>.mp4` |
| the site | `publish.py <film> <poster time>` | `docs/_static/video/` |

The picture is a web page (`film/`) that draws any frame on request (`window.FILM.frame(t)`);
`render.py` steps it in headless Chromium and encodes the frames with the soundtrack. `voice.py` runs
in an environment with `chatterbox-tts` and `torchaudio`; `render.py` needs `playwright`.
