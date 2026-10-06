/*
 * MetaPhlAn profile helpers (ADR-0018)
 *
 * Profiles are defined in conf/metaphlan_profiles.config; this is the one place
 * the "does HUMAnN reuse the main pass?" rule lives, so DATABASES, the HUMANN
 * subworkflow and the manifest cannot disagree about it.
 */

// True when the selected HUMAnN bundle names the same MetaPhlAn profile as the
// main taxonomy pass: HUMAnN then consumes the main full-depth profile and
// neither a second MetaPhlAn index install nor a second MetaPhlAn run happens.
def humann_reuses_main_metaphlan() {
    return params.humann_bundles[params.humann_bundle].metaphlan_profile == params.metaphlan_profile
}
