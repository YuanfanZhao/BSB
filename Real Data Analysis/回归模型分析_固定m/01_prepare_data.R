# ============================================================================
# 01_prepare_data.R --- 准备"固定 m + 多协变量"的优质真实数据集（回归分析用）
# ----------------------------------------------------------------------------
# 输出 data/*.csv（列：x = 成功次数, m = 固定试验次数, 协变量...）
#   1) MathExam14W（Innsbruck 大学数学101期末考试，psychotools 包）
#      729 名学生，13 道单选题，x = 答对题数（0~13），m = 13 固定。
#      协变量（5 个）：tests、gender、study、semester、attempt。
#   2) PISA 2000 Reading（CDM 包 data.pisa00R.ct，Chen & de la Torre 2014）
#      1095 名学生，23 道二分阅读题，x = 答对题数（0~23），m = 23 固定。
#      协变量（3 个）：female、HISEI（家庭社会经济指数）、AGE（月龄）。
#   3) CTB（CTB/McGraw-Hill 学业成就测验，structree 包）——探索用
#      1500 名学生，56 道多选题，x = 答对题数（21~46），m = 56 固定。
#      协变量（5 个）：gender、type、size、bachelor、language。
#      （注：CTB 条件分布接近二项、无明显过度离散，Bern-Bino 在其上不占优，
#       仅作对照探索；正式实例以 1)、2) 为准。）
# 使用 .libPaths() 指向用户库（系统库 D:/R-4.5.2/library 不可写）。
# ============================================================================
.libPaths(c("C:/Users/Zhao Yuanfan/AppData/Local/R/win-library/4.5", .libPaths()))
dir.create("data", showWarnings = FALSE)
suppressMessages({library(psychotools); library(CDM); library(structree)})

# ---- 1. MathExam14W ----
data("MathExam14W", package = "psychotools")
d1 <- data.frame(
  x        = MathExam14W$nsolved,        # 答对题数（0~13）
  m        = 13L,                        # 固定试验次数
  tests    = MathExam14W$tests,          # 考前在线练习完成数（9~26）
  gender   = factor(MathExam14W$gender), # female / male
  study    = factor(MathExam14W$study),  # 155 = 4 年制, 571 = 3 年制
  semester = MathExam14W$semester,       # 已读学期数（1~21）
  attempt  = factor(MathExam14W$attempt) # 考试尝试次数（1~5）
)
write.csv(d1, "data/mathexam_data.csv", row.names = FALSE)
cat("MathExam14W: n =", nrow(d1), " m =", unique(d1$m),
    " x range =", paste(range(d1$x), collapse = "-"), "\n")

# ---- 2. PISA 2000 Reading ----
e <- new.env(); data("data.pisa00R.ct", package = "CDM", envir = e)
dp <- e$data.pisa00R.ct$data
itemnms <- grep("^R", names(dp), value = TRUE)
bin <- vapply(itemnms, function(nm) {
  z <- dp[[nm]]; is.numeric(z) && all(z %in% c(0, 1), na.rm = TRUE) && !any(is.na(z))
}, logical(1))
items <- itemnms[bin]
d2 <- data.frame(
  x      = rowSums(as.matrix(dp[, items])), # 答对题数（0~23）
  m      = length(items),                   # 固定试验次数（23）
  female = dp$female,                       # 0/1
  hisei  = dp$HISEI,                        # 家庭社会经济指数（16~99）
  age    = dp$AGE                           # 月龄
)
write.csv(d2, "data/pisa00_data.csv", row.names = FALSE)
cat("PISA2000Read: n =", nrow(d2), " m =", unique(d2$m),
    " x range =", paste(range(d2$x), collapse = "-"), "\n")

# ---- 3. CTB（探索用对照） ----
data("CTB", package = "structree")
d3 <- data.frame(
  x        = CTB$score,
  m        = 56L,
  gender   = CTB$gender,
  type     = factor(CTB$type),
  size     = CTB$size,
  bachelor = CTB$bachelor,
  language = CTB$language
)
write.csv(d3, "data/ctb_data.csv", row.names = FALSE)
cat("CTB: n =", nrow(d3), " m =", unique(d3$m),
    " x range =", paste(range(d3$x), collapse = "-"), "\n")
cat("data files written.\n")
