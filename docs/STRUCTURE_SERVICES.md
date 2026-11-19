# Structure-service access for a small STOP4 trial

Checked 2026-11-19. Prices and access rules are provider-controlled; the prototype
does not create accounts, buy credits, deploy endpoints or submit predictions.

| Service | Actual method | Access and practical fit |
| --- | --- | --- |
| AlphaFold Server | Google DeepMind AlphaFold 3 | Free non-commercial web service with Google login. EMBL-EBI documents 30 jobs/day. Good first route for a few true AF3 predictions; a supported public submission API was not confirmed. |
| Google Cloud AF3 | Google DeepMind AlphaFold 3 | Official prediction API, including end-to-end sequence/MSA/inference. Current documentation requires allowlist/commercial subscription access, an endpoint and output storage. Greater setup than a first academic trial; no simple per-protein price verified. |
| Neurosnap | AlphaFold2, Boltz-2, Protenix, Chai-1 and other named tools | Documented API with an account API key and job/file operations. Pay-as-you-go lists $1/credit, $0.90 for verified academics, including full API access. Charges depend on tool/runtime, not one fixed protein price. The free plan includes one introductory credit, but free-tier API/model eligibility needs account confirmation. Its “AlphaFold3” catalogue labels include reproductions; preserve actual model names. |
| Tamarind Bio | AlphaFold2/ColabFold, Boltz, Chai and other tools | Documented API; free page lists 10 jobs/month, academic extensions are available. Open API access is explicitly listed under Premium; academic API entitlement and price need confirmation. The provider's “AlphaFold3 commercially” article discusses Boltz/Chai alternatives. |
| Foldify / Foldify-Open | Includes actual AlphaFold 3 | The published open service is a possible free academic trial route. The institutional manual describes e-INFRA login/group access, while the open repository supports guest sessions. Public submission-API support and current open-service availability were not verified. |

For a few original AF3 structures, start with AlphaFold Server and the generated
`.server.json` uploads. For an API-first prototype without managing GPU/cloud
infrastructure, Neurosnap pay-as-you-go is a candidate **if a named alternative
such as Boltz-2 is acceptable**. Run one small, explicitly selected job and inspect
actual credit usage before estimating a whole multi-species study. Google Cloud
is the confirmed official AF3 API route, with its access/deployment requirements.

## Verified primary sources

- [Google DeepMind AF3 source and access routes](https://github.com/google-deepmind/alphafold3)
- [Official Google Cloud AF3 endpoint guide](https://docs.cloud.google.com/gemini-enterprise-agent-platform/models/open-models/alphafold-3)
- [EMBL-EBI: AlphaFold Server access, free use and quota](https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/alphafold-server-your-gateway-to-alphafold-3/)
- [Neurosnap catalogue](https://neurosnap.ai/services), [API](https://neurosnap.ai/docs/automation/api-reference), [pricing and credit definitions](https://neurosnap.ai/)
- [Tamarind pricing](https://www.tamarind.bio/pricing), [API](https://app.tamarind.bio/api-docs/overview), [its AF3/alternative distinction](https://www.tamarind.bio/blog/how-to-use-alphafold3-commercially)
- [Foldify institutional manual](https://docs.cerit.io/en/docs/web-apps/foldify), [open-service source](https://github.com/CERIT-SC/foldify-open)
- [AFDB API](https://alphafold.ebi.ac.uk/api-docs), [current response-field changes](https://www.ebi.ac.uk/pdbe/news/breaking-changes-afdb-predictions-api)
