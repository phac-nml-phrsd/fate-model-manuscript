# Ottawa neighbourhood–WWTP benchmark data

## Provenance and scope

`Ottawa_data.xlsx` contains aggregate wastewater measurements supplied by
the Ottawa field-study collaborators (University of Carleton, Ottawa, Ontario, Canada)
and used as an independent benchmark in
the revised fate-model manuscript. The
manuscript credits Lena Carolin Bitter and Banu Ormecy with providing the Ottawa field
data and states that data underlying this comparison are publicly available
in this repository.

The sampling campaign compared an upper-sewershed neighbourhood with the
Robert O. Pickard Environmental Centre. The paired calculation uses
2021-03-22 through 2021-05-18. The manuscript describes the principal paired
measurements as spanning April 8 through May 27, 2021; the workbook contains
additional surrounding context dates. Dates used in calculations are retained
in the published paired-ratio table.

The workbook contains aggregate environmental measurements only. No names,
addresses, clinical records, or individual-level health information were
identified during this repository review.

## Redistribution status

The revised manuscript designates these data for public release, and the
workbook is held in the manuscript authors' working repository. The workbook
itself does not contain a licence or written data-owner authorization. Thus,
the repository supports intended public release, but final confirmation of
redistribution authority must be made by the project lead/data owners before
publication.

## Workbook sheets and fields

Sheets are `N1_UF`, `N1_CF`, `N2_UF`, `N2_CF`, `PMMoV_UF`, `PMMoV_CF`, `TS`,
`VS`, `flow`, and `Residence Time`.

- `date`: sampling date.
- `WWTP`, `D`: treatment-plant and neighbourhood-D measurements used in the
  revised benchmark. 
- `target`: measurement and sample-fraction identifier.
- `note`: sampling note, including the documented grab sample.
- N1 and N2: SARS-CoV-2 nucleocapsid targets, gene copies per mL.
- PMMoV: pepper mild mottle virus, gene copies per mL.
- TS: total solids, mg/L.
- VS: volatile solids, mg/L.
- flow: wastewater flow rate, L/s.
- HRT: hydraulic residence time, hours.

Fractions are `UF` (liquid, ultrafiltration), `CF` (solids, centrifugation),
and derived `TOTAL` (paired UF + CF concentrations).

## Processing and exclusions

The analysis reads every measurement sheet except sheets containing
`residence`, `HRT`, `travel`, or `retention`; HRT is contextual rather than an
apparent-loss input. Duplicate date/site/target/fraction values are reduced to
their median. Ratios require positive, non-missing paired WWTP and site-D
values. Zero and missing values therefore do not enter log10 ratios or
normalization denominators. The public workbook and derived outputs contain
only the WWTP and neighbourhood-D analysis.

The grab-sample observation remains flagged in `note` and is not automatically
excluded. This benchmark is contextual and was not used for model calibration
or direct validation.
