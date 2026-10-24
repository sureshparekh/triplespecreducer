# The faces BRAIN draws its interface in

Vendored rather than assumed. A font that is not on a machine falls back
silently, and an interface whose type silently changes is one nobody designed.
These are loaded from this folder at startup, so BRAIN looks the same on a
machine that has never seen them. They are the faces of the documentation
website, so the interface and the website read as one thing.

| File | Used for | Licence |
|---|---|---|
| `HankenGrotesk-*.ttf` | the interface: labels, buttons, fields, prose | SIL OFL 1.1 |
| `SourceSerif4-*.ttf` | the wordmark and the title of each step | SIL OFL 1.1 |
| `JetBrainsMono-*.ttf` | every number, with tabular figures | SIL OFL 1.1 |

The full licence text for each is in the matching `OFL-*.txt`.

Each file is one weight cut from the family's variable font in the
google/fonts repository (Regular, Medium, SemiBold and Bold for Hanken
Grotesk; Light, Regular and SemiBold for Source Serif 4; Regular and Medium
for JetBrains Mono), so that Qt, which handles static weights more reliably
than variable axes, finds every weight the stylesheet asks for.
