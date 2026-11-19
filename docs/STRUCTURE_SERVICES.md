# Structure-service access for a small STOP4 trial

Checked 2026-11-19. Prices and access rules are provider-controlled; the prototype
does not create accounts, buy credits, deploy endpoints or submit predictions.

| Service | Actual method | Access and practical fit |
| --- | --- | --- |
| AlphaFold Server | Google DeepMind AlphaFold 3 | Free non-commercial web service with Google login. EMBL-EBI documents 30 jobs/day. Good first route for a few true AF3 predictions; a supported public submission API was not confirmed. |
| Google Cloud AF3 | Google DeepMind AlphaFold 3 | Official prediction API, including end-to-end sequence/MSA/inference. Current documentation requires allowlist/commercial subscription access, an endpoint and output storage. Greater setup than a first academic trial; no simple per-protein price verified. |
| NVIDIA-hosted OpenFold3 NIM | OpenFold Consortium OpenFold3, a reimplementation of the AF3 approach | A hosted prediction endpoint and Get API Key flow are documented. NVIDIA's FAQ offers free NIM API access for Developer Program prototyping. Register, generate a key, then check account/model eligibility and current rate limits; no guaranteed number of protein jobs was verified. Suitable first API trial if an AF3 reimplementation is acceptable. Hosted trial pages state that input/output can be recorded; use public controls first. |
| Neurosnap | AlphaFold2, Boltz-2, Protenix, Chai-1 and other named tools | Documented API with an account API key and job/file operations. Pay-as-you-go lists $1/credit, $0.90 for verified academics, including full API access. Charges depend on tool/runtime, not one fixed protein price. The free plan includes one introductory credit, but free-tier API/model eligibility needs account confirmation. Its “AlphaFold3” catalogue labels include reproductions; preserve actual model names. |
| Tamarind Bio | AlphaFold2/ColabFold, Boltz, Chai and other tools | Documented API; free page lists 10 jobs/month, academic extensions are available. Open API access is explicitly listed under Premium; academic API entitlement and price need confirmation. The provider's “AlphaFold3 commercially” article discusses Boltz/Chai alternatives. |
| Foldify / Foldify-Open | Includes actual AlphaFold 3 | The published open service is a possible free academic trial route. The institutional manual describes e-INFRA login/group access, while the open repository supports guest sessions. Public submission-API support and current open-service availability were not verified. |

The project's preferred route for new predictions is **AlphaFold Server**. Use
the generated `upload.server.json` batch, run its drafts, then retain and import
each downloaded result ZIP with `af3-import`. A real trial on the saved POU5F3
proteins of Gallus gallus and Pavo cristatus completed on 2026-11-19: two jobs,
one seed each, five returned samples per job. These are new Server AF3 predictions,
separate from earlier ready-made AFDB controls. These experimental jobs do not
approve either gene model as a final reference. After the user's clarification,
further exports require the separate **«Утвердить модель как эталон»** action in
the STOP3 report; a reviewed choice alone does not grant approval;
candidate selections and alignments cannot authorize submission. See `docs/STAGE4.md` for the
request/result boundary and later provider adapters.

For an API-first prototype without managing GPU/cloud
infrastructure, first check NVIDIA's hosted OpenFold3 trial **if a reimplementation
is acceptable**. Its hosted endpoint does not require a local Docker/GPU setup.
Neurosnap pay-as-you-go is another candidate for named alternatives such as
Boltz-2. Run one small, explicitly selected job and inspect actual limits/usage
before estimating a whole multi-species study. Google Cloud is the confirmed
official AF3 API route, with its access/deployment requirements.

## First hosted API access

1. Open the [NVIDIA OpenFold3 page](https://build.nvidia.com/openfold/openfold3)
   and choose **Get API Key**; register/sign in to the NVIDIA Developer Program.
2. Confirm the account can use OpenFold3 and inspect current trial limits.
   The documented endpoint is
   `https://health.api.nvidia.com/v1/biology/openfold/openfold3/predict`.
3. Keep the key in a local environment/credential store, outside Git and chat.
   Select a public reference sequence for the first authorized smoke test.
   This project has not submitted such a test or configured this provider yet.

## Verified primary sources

- [Google DeepMind AF3 source and access routes](https://github.com/google-deepmind/alphafold3)
- [Official Google Cloud AF3 endpoint guide](https://docs.cloud.google.com/gemini-enterprise-agent-platform/models/open-models/alphafold-3)
- [EMBL-EBI: AlphaFold Server access, free use and quota](https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/alphafold-server-your-gateway-to-alphafold-3/)
- [NVIDIA OpenFold3 hosted trial](https://build.nvidia.com/openfold/openfold3), [model/API reference](https://docs.api.nvidia.com/nim/reference/openfold-openfold3), [Developer Program API access FAQ](https://docs.api.nvidia.com/nim/docs/product)
- [Neurosnap catalogue](https://neurosnap.ai/services), [API](https://neurosnap.ai/docs/automation/api-reference), [pricing and credit definitions](https://neurosnap.ai/)
- [Tamarind pricing](https://www.tamarind.bio/pricing), [API](https://app.tamarind.bio/api-docs/overview), [its AF3/alternative distinction](https://www.tamarind.bio/blog/how-to-use-alphafold3-commercially)
- [Foldify institutional manual](https://docs.cerit.io/en/docs/web-apps/foldify), [open-service source](https://github.com/CERIT-SC/foldify-open)
- [AFDB API](https://alphafold.ebi.ac.uk/api-docs), [current response-field changes](https://www.ebi.ac.uk/pdbe/news/breaking-changes-afdb-predictions-api)
