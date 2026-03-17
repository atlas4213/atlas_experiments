#!/usr/bin/env python3
"""
01_engagement_analysis.py
--------------------------
Core analysis for the Healthcare Product Usage & Engagement project.

Analyzes:
  1. Funnel analysis  — onboarding → active → power user conversion
  2. Feature adoption — which features drive engagement and retention
  3. Churn analysis   — leading indicators of disengagement
  4. Cohort analysis  — retention curves by join month
  5. Outcome correlation — engagement tier vs clinical outcomes
  6. Department benchmarking — cross-department usage patterns

Outputs results/ JSON and TSV files consumed by the dashboard and report.
"""

import csv, json, math, statistics
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

DATA = Path("data")
OUT  = Path("results")
OUT.mkdir(exist_ok=True)

# ── Loaders ────────────────────────────────────────────────────────────────────
def load_csv(fname):
    with open(DATA / fname) as f:
        return list(csv.DictReader(f))

print("Loading data...")
users    = load_csv("users.csv")
sessions = load_csv("sessions.csv")
events   = load_csv("feature_events.csv")
outcomes = load_csv("outcomes.csv")
print(f"  {len(users)} users | {len(sessions)} sessions | {len(events)} events")

user_map    = {u["user_id"]: u for u in users}
outcome_map = {o["user_id"]: o for o in outcomes}

# ── 1. FUNNEL ANALYSIS ─────────────────────────────────────────────────────────
print("\n[1] Funnel Analysis")

total_invited    = 600   # simulated invites sent
total_registered = len(users)

# Activated = had at least 1 session in first 7 days
session_by_user = defaultdict(list)
for s in sessions:
    session_by_user[s["user_id"]].append(s["session_date"])

activated = 0
for user in users:
    join = datetime.strptime(user["join_date"], "%Y-%m-%d")
    user_sessions = sorted(session_by_user[user["user_id"]])
    if user_sessions:
        first = datetime.strptime(user_sessions[0], "%Y-%m-%d")
        if (first - join).days <= 7:
            activated += 1

active_users    = sum(1 for u in users if u["engagement_tier"] in ["Active", "Power User"])
power_users     = sum(1 for u in users if u["engagement_tier"] == "Power User")
churned_users   = sum(1 for u in users if u["churned"] == "Yes")

funnel = {
    "stages": [
        {"stage": "Invited",       "count": total_invited,    "pct_of_top": 100.0},
        {"stage": "Registered",    "count": total_registered, "pct_of_top": round(total_registered/total_invited*100,1)},
        {"stage": "Activated",     "count": activated,        "pct_of_top": round(activated/total_invited*100,1)},
        {"stage": "Active User",   "count": active_users,     "pct_of_top": round(active_users/total_invited*100,1)},
        {"stage": "Power User",    "count": power_users,      "pct_of_top": round(power_users/total_invited*100,1)},
    ],
    "drop_off": []
}
for i in range(1, len(funnel["stages"])):
    prev = funnel["stages"][i-1]["count"]
    curr = funnel["stages"][i]["count"]
    funnel["drop_off"].append({
        "from": funnel["stages"][i-1]["stage"],
        "to":   funnel["stages"][i]["stage"],
        "lost": prev - curr,
        "conversion_rate": round(curr/prev*100, 1)
    })

print(f"  Funnel: {total_invited} invited → {total_registered} registered → "
      f"{activated} activated → {active_users} active → {power_users} power users")

# ── 2. FEATURE ADOPTION ────────────────────────────────────────────────────────
print("\n[2] Feature Adoption")

feature_users   = defaultdict(set)
feature_time    = defaultdict(list)
feature_complete= defaultdict(int)
feature_total   = defaultdict(int)

for e in events:
    feature_users[e["feature"]].add(e["user_id"])
    feature_time[e["feature"]].append(int(e["time_spent_s"]))
    feature_total[e["feature"]] += 1
    if e["completed"] == "Yes":
        feature_complete[e["feature"]] += 1

feature_adoption = []
for feat in sorted(feature_users, key=lambda x: -len(feature_users[x])):
    times = feature_time[feat]
    feature_adoption.append({
        "feature":          feat,
        "unique_users":     len(feature_users[feat]),
        "adoption_rate":    round(len(feature_users[feat]) / len(users) * 100, 1),
        "total_events":     feature_total[feat],
        "avg_time_s":       round(statistics.mean(times)),
        "median_time_s":    round(statistics.median(times)),
        "completion_rate":  round(feature_complete[feat] / feature_total[feat] * 100, 1),
    })
    print(f"  {feat:<30} {len(feature_users[feat]):>4} users "
          f"({round(len(feature_users[feat])/len(users)*100,1)}%)")

# ── 3. CHURN ANALYSIS ─────────────────────────────────────────────────────────
print("\n[3] Churn Analysis")

churn_by_tier = defaultdict(lambda: {"churned":0,"total":0})
for u in users:
    t = u["engagement_tier"]
    churn_by_tier[t]["total"] += 1
    if u["churned"] == "Yes":
        churn_by_tier[t]["churned"] += 1

churn_analysis = {
    "overall_churn_rate": round(churned_users / len(users) * 100, 1),
    "by_tier": [
        {
            "tier": t,
            "total": v["total"],
            "churned": v["churned"],
            "churn_rate": round(v["churned"]/v["total"]*100, 1)
        }
        for t, v in sorted(churn_by_tier.items(),
                           key=lambda x: -x[1]["churned"]/x[1]["total"])
    ]
}

churn_by_dept = defaultdict(lambda: {"churned":0,"total":0})
for u in users:
    churn_by_dept[u["department"]]["total"] += 1
    if u["churned"] == "Yes":
        churn_by_dept[u["department"]]["churned"] += 1

churn_analysis["by_department"] = [
    {"dept": d, "churn_rate": round(v["churned"]/v["total"]*100,1), "total": v["total"]}
    for d, v in sorted(churn_by_dept.items(), key=lambda x: -x[1]["churned"]/x[1]["total"])
]

print(f"  Overall churn: {churn_analysis['overall_churn_rate']}%")
for row in churn_analysis["by_tier"]:
    print(f"  {row['tier']:<16}: {row['churn_rate']}% churn rate")

# ── 4. COHORT RETENTION ────────────────────────────────────────────────────────
print("\n[4] Cohort Retention")

cohorts = defaultdict(list)
for u in users:
    month = u["join_date"][:7]
    cohorts[month].append(u["user_id"])

monthly_active = defaultdict(set)
for s in sessions:
    month = s["session_date"][:7]
    monthly_active[month].add(s["user_id"])

cohort_retention = []
all_months = sorted(cohorts.keys())[:6]  # Jan–Jun cohorts
for c_month in all_months:
    c_users = set(cohorts[c_month])
    c_idx   = all_months.index(c_month)
    row = {"cohort": c_month, "size": len(c_users), "retention": []}
    for offset in range(7):
        future_months = sorted(monthly_active.keys())
        target_months = [m for m in future_months if m >= c_month]
        if offset < len(target_months):
            m = target_months[offset]
            active = len(c_users & monthly_active[m])
            row["retention"].append(round(active/len(c_users)*100, 1))
        else:
            row["retention"].append(None)
    cohort_retention.append(row)
    print(f"  Cohort {c_month}: {len(c_users)} users, M0={row['retention'][0]}%")

# ── 5. OUTCOME CORRELATIONS ────────────────────────────────────────────────────
print("\n[5] Outcome Correlations by Engagement Tier")

outcome_by_tier = defaultdict(lambda: defaultdict(list))
for o in outcomes:
    uid  = o["user_id"]
    tier = user_map[uid]["engagement_tier"]
    for k in ["care_gap_closure_rate", "avg_alert_response_min",
              "patient_satisfaction_score", "documentation_accuracy_pct",
              "30day_readmission_rate"]:
        outcome_by_tier[tier][k].append(float(o[k]))

outcome_summary = []
tier_order = [t for t in ["Power User", "Active", "Casual", "At-Risk"] if outcome_by_tier[t]["care_gap_closure_rate"]]
for tier in tier_order:
    data = outcome_by_tier[tier]
    outcome_summary.append({
        "tier": tier,
        "n": len(data["care_gap_closure_rate"]),
        "care_gap_closure_pct":     round(statistics.mean(data["care_gap_closure_rate"])*100, 1),
        "avg_alert_response_min":   round(statistics.mean(data["avg_alert_response_min"]), 1),
        "patient_satisfaction":     round(statistics.mean(data["patient_satisfaction_score"]), 2),
        "documentation_accuracy":   round(statistics.mean(data["documentation_accuracy_pct"]), 1),
        "readmission_rate_pct":     round(statistics.mean(data["30day_readmission_rate"])*100, 1),
    })
    print(f"  {tier:<16}: satisfaction={outcome_summary[-1]['patient_satisfaction']}, "
          f"readmission={outcome_summary[-1]['readmission_rate_pct']}%, "
          f"care gap closure={outcome_summary[-1]['care_gap_closure_pct']}%")

# ── 6. DEPARTMENT BENCHMARKING ─────────────────────────────────────────────────
print("\n[6] Department Benchmarking")

dept_stats = defaultdict(lambda: {"sessions":0,"users":set(),"events":0,"features":set()})
for s in sessions:
    dept = user_map[s["user_id"]]["department"]
    dept_stats[dept]["sessions"] += 1
    dept_stats[dept]["users"].add(s["user_id"])

for e in events:
    dept = user_map[e["user_id"]]["department"]
    dept_stats[dept]["events"] += 1
    dept_stats[dept]["features"].add(e["feature"])

dept_benchmark = []
for dept, v in sorted(dept_stats.items()):
    n_users = sum(1 for u in users if u["department"] == dept)
    power   = sum(1 for u in users if u["department"] == dept and u["engagement_tier"] == "Power User")
    dept_benchmark.append({
        "department":       dept,
        "total_users":      n_users,
        "active_users":     len(v["users"]),
        "power_user_pct":   round(power/n_users*100, 1),
        "sessions_per_user":round(v["sessions"]/max(1,len(v["users"])), 1),
        "features_adopted": len(v["features"]),
        "avg_events_per_session": round(v["events"]/max(1,v["sessions"]), 1),
    })
    print(f"  {dept:<16}: {len(v['users'])} active users, "
          f"{round(v['sessions']/max(1,len(v['users'])),1)} sessions/user")

# ── 7. MONTHLY TREND ──────────────────────────────────────────────────────────
monthly_sessions = defaultdict(int)
monthly_users    = defaultdict(set)
for s in sessions:
    m = s["session_date"][:7]
    monthly_sessions[m] += 1
    monthly_users[m].add(s["user_id"])

monthly_trend = [
    {"month": m, "sessions": monthly_sessions[m], "active_users": len(monthly_users[m])}
    for m in sorted(monthly_sessions.keys())
]

# ── Save all results ───────────────────────────────────────────────────────────
results = {
    "funnel":            funnel,
    "feature_adoption":  feature_adoption,
    "churn":             churn_analysis,
    "cohort_retention":  cohort_retention,
    "outcomes":          outcome_summary,
    "dept_benchmark":    dept_benchmark,
    "monthly_trend":     monthly_trend,
    "summary_stats": {
        "total_users":    len(users),
        "total_sessions": len(sessions),
        "total_events":   len(events),
        "avg_session_min": round(statistics.mean(float(s["duration_min"]) for s in sessions), 1),
        "power_user_pct": round(power_users/len(users)*100, 1),
        "churn_rate_pct": round(churned_users/len(users)*100, 1),
    }
}

with open(OUT / "analysis_results.json", "w") as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to results/analysis_results.json")
print("Analysis complete.")
